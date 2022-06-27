from pandas import read_csv, DataFrame
from os.path import join
import numpy as np
import json
from tqdm import tqdm

from cobra.io import write_sbml_model

from troppo.methods.reconstruction.fastcore import FastcoreProperties
from troppo.methods_wrappers import ReconstructionWrapper


class ReconstructModel(object):

    def __init__(self, generic_model, raw_counts_file, geneid_mapping_file, raw_counts_cols=None,
                 tas_global_min=0.1, tas_local=.25, tas_global_max=.75, flux_threshold=1e-7, solver='CPLEX'):
        # Model that will be used for reconstruction:
        self.generic_model = generic_model.copy()
        # Mapping between gene symbols and ensemle ids:
        with open(geneid_mapping_file) as f: self.geneID_mapping = json.load(f)
        # Tas thresholds:
        self.tas_thresholds = {'global_min': tas_global_min, 'local': tas_local, 'global_max': tas_global_max}
        # Raw, processed and tas expression counts are stored in dictionary:
        rc = read_csv(raw_counts_file, header=0, index_col=0)
        self.expression = {}
        if raw_counts_cols is None:
            self.expression = {'raw_counts': rc}
        else:
            common = list(set(raw_counts_cols) & set(rc.columns.to_list()))
            if len(common) == 0:
                ValueError('Invalid raw_counts_cols.')
            else:
                self.expression = {'raw_counts': rc[common]}
        self._cpm()
        # TAS, and RAS, scores for only the genes in the generic model will be stored in a dictionary, as well as
        # FASTCORE input parameters will be stored here:
        self.reconstruction_inputs = {}
        self._set_fastcore_parameters(flux_threshold, solver)
        # Name of the future reconstructed model:
        self.reconstructed_models_names = self.expression['raw_counts'].columns
        # Where fastcore properties will be stored:
        self.fastcore_properties = {}
        # Where the reconstructed model will be saved:
        self.reconstructed_models = {}

    def _cpm(self):
        self.expression['processed_counts'] = self.expression['raw_counts'] / self.expression['raw_counts'].sum(0) * 1e6

    def _calculate_tas(self):
        global_min, global_max = np.quantile(self.expression['processed_counts'], [.1, .75])
        local_thresholds = np.quantile(self.expression['processed_counts'], [.25], axis=1)[0]

        tas_counts = np.array(self.expression['processed_counts'].copy())
        active = tas_counts >= global_max
        inactive = tas_counts <= global_min
        moderate = np.where(np.logical_and(tas_counts > global_min, tas_counts < global_max))
        tas_counts[active] = np.log2(tas_counts[active] / global_max)
        tas_counts[inactive] = np.log2(tas_counts[inactive] / global_min)
        for gene, col in zip(moderate[0], moderate[1]):
            if local_thresholds[gene] == 0:
                tas_counts[gene, col] = 1
            else:
                tas_counts[gene, col] = np.log2(tas_counts[gene, col] / local_thresholds[gene])

        self.expression['tas_counts'] = DataFrame(tas_counts.copy(),
                                                  index=self.expression['raw_counts'].index,
                                                  columns=self.expression['raw_counts'].columns)

    def _convert_gene_ids(self):
        model_ids = [i.id for i in self.generic_model.genes]
        from_ids = self.expression['tas_counts'].index
        counts = self.expression['tas_counts'].copy()
        to_ids = np.array([self.geneID_mapping.get(key) for key in from_ids])
        counts = counts.loc[[x is not None for x in to_ids],]
        counts.index = [x[0] for x in list(filter(None, to_ids))]
        no_info_genes = np.setdiff1d(model_ids, counts.index)
        counts = counts.append(DataFrame(0, index=no_info_genes, columns=counts.columns))
        self.reconstruction_inputs['TAS'] = counts.copy()

    def _get_ras(self):
        ras = np.zeros((len(self.generic_model.reactions), self.reconstruction_inputs['TAS'].shape[1]))
        rxn_names = []
        cols = self.reconstruction_inputs['TAS'].columns
        for i, reaction in tqdm(enumerate(self.generic_model.reactions)):
            gpr_reaction = reaction.gene_reaction_rule
            rxn_names = rxn_names + [reaction.id]
            if gpr_reaction is not '':
                or_vals = {i: [] for i in cols}
                for or_part in gpr_reaction.strip().replace('(', '').replace(')', '').split(' or '):
                    and_part = [x.strip() for x in or_part.strip().split(' and ')]
                    for col in cols:
                        or_vals[col] = or_vals[col] + [self.reconstruction_inputs['TAS'].loc[and_part, col].min()]
                for j, col in enumerate(cols):
                    ras[i, j] = np.array(or_vals[col]).max()
        self.reconstruction_inputs['RAS'] = DataFrame(ras.copy(), index=rxn_names, columns=cols)

    def _get_core_reactions(self, media):
        core_reactions_id = {}
        core_reactions_idx = {}
        for col in self.reconstruction_inputs['RAS'].columns:
            # Get ID of core reactions:
            core_reactions_id[col] = [k for (k, v) in self.reconstruction_inputs['RAS'][col].items() if v > 0]
            core_reactions_id[col] = np.unique(core_reactions_id[col] + media).tolist()
            # Get the indexes of said reactions:
            core_reactions_idx[col] = [self.generic_model.reactions.index(rxn) for rxn in core_reactions_id[col]]
            core_reactions_idx[col] = np.unique(core_reactions_idx[col] + [self.generic_model.reactions.index(rxn)
                                                                           for rxn in media]).tolist()
        self.reconstruction_inputs['core_reactions_id'] = core_reactions_id
        self.reconstruction_inputs['core_reactions_idx'] = core_reactions_idx

    def _set_fastcore_parameters(self, flux_threshold, solver):
        self.reconstruction_inputs['flux_threshold'] = flux_threshold
        self.reconstruction_inputs['solver'] = solver

    def _set_algorithm_properties(self):
        for col in self.reconstructed_models_names:
            self.fastcore_properties[col] = FastcoreProperties(self.reconstruction_inputs['core_reactions_idx'][col],
                                                               flux_threshold=
                                                               self.reconstruction_inputs['flux_threshold'],
                                                               solver=self.reconstruction_inputs['solver'])

    def run(self, media):
        # Calculate TAS:
        print('Calculating TAS scores...')
        self._calculate_tas()
        self._convert_gene_ids()
        # Calculate RAS:
        print('\n\nCalculating RAS scores...')
        self._get_ras()
        print('\n\nGetting core reactions from RAS scores...')
        self._get_core_reactions(media)
        # Set properties for each reconstruction:
        print('\n\nSetting properties for reconstruction(s)...')
        self._set_algorithm_properties()
        # Reconstruct models for each column in raw_counts file:
        print('\n\nReconstructing model(s): ')
        for col in self.reconstructed_models_names:
            print('\n--->', col)
            x = ReconstructionWrapper(self.generic_model)
            res = x.run(self.fastcore_properties[col])
            new_model = self.generic_model.copy()
            for r_idx in range(0, len(self.generic_model.reactions)):
                if r_idx not in res:
                    new_model.reactions[r_idx].knock_out()
            self.reconstructed_models[col] = new_model
        print('Done!')

    def retrieve_data_info(self):
        if 'tas_counts' not in self.expression.keys():
            self._calculate_tas()
            self._convert_gene_ids()
            self._get_ras()
        data_info = {'CPM': self.expression['processed_counts'],
                     'TAS': self.reconstruction_inputs['TAS'],
                     'RAS': self.reconstruction_inputs['RAS']}
        return data_info


