from pandas import read_csv
import warnings

from troppo.methods.reconstruction.fastcore import FastcoreProperties
from troppo.methods.reconstruction.gimme import GIMMEProperties
from troppo.methods.reconstruction.imat import IMATProperties
from troppo.methods.reconstruction.tINIT import tINITProperties
from troppo.methods.reconstruction.corda import CORDAProperties

from troppo.methods_wrappers import ModelBasedWrapper, ReconstructionWrapper


class ReconstructModel(object):

    def __init__(self, generic_model, algorithm, reaction_calls_file, reconstructed_model_name, **kwargs):
        # Model that will be useobjectivesd for reconstruction:
        self.generic_model = generic_model.copy()
        # Check if algorithm wanted is available:
        if not all(elem in ['Fastcore', 'GIMME', 'IMAT', 'tINIT', 'CORDA'] for elem in [algorithm]):
            warnings.warn('algorithm given not available. Algorithms available: Fastcore, GIMME, IMAT, tINIT, CORDA')
            return
        else:
            self.algorithm = algorithm
        # Get the reactions calls for the wanted algorithm:
        self.reaction_calls = self.get_reaction_calls_from_file(reaction_calls_file)
        # Name of the future reconstructed model:
        self.reconstructed_model_name = reconstructed_model_name
        # Get the algorithms' properties:
        self.algorithm_properties = self.get_algorithm_properties(**kwargs)
        # Where the recnstructed model will be saved:
        self.reconstructed_model = None

    def get_reaction_calls_from_file(self, reaction_calls_file):
        return read_csv(reaction_calls_file, header=0, index_col=0).to_dict()[self.algorithm]

    def get_algorithm_properties(self, **kwargs):

        # Check if all necessary arguments for the algorithm wanted are given. If so, create the properties object:
        if self.algorithm == 'FASTcore':
            if not all(elem in kwargs.keys() for elem in ['flux_threshold', 'solver']):
                warnings.warn('FASTcore properties missing. The following arguments must be present: '
                              'flux_threshold, solver')
                return None
            else:
                # Get 'Present' reactions (id):
                core_reactions_id = [k for (k, v) in self.reaction_calls.items() if v == 'Present']
                # Get the indexes of said reactions:
                core_reactions_idx = [self.generic_model.reactions.index(rxn) for rxn in core_reactions_id]
                # Create FASTcoreProperties object:
                return FastcoreProperties(core_reactions_idx, flux_threshold=kwargs.get('flux_threshold'),
                                          solver=kwargs.get('solver'))

        elif self.algorithm == 'GIMME':
            if not all(elem in kwargs.keys() for elem in ['objectives', 'obj_frac', 'flux_threshold']):
                warnings.warn('GIMME properties missing. The following arguments must be present: objectives, '
                              'obj_frac, flux_threshold')
                return None
            else:
                # Convert list of objectives with rereconstructed_modelaction ids to reaction indexes:
                objectives = kwargs.get('objectives')
                for x_idx, x_val in enumerate(objectives):
                    for y_idx, y_val in enumerate(x_val):
                        objectives[x_idx][y_idx] = self.generic_model.reactions.index(y_val)
                # Create GIMMEProperties object:
                return GIMMEProperties(exp_vector=list(self.reaction_calls.values()), objectives=objectives,
                                       obj_frac=kwargs.get('obj_frac'), preprocess=False,
                                       flux_threshold=kwargs.get('flux_threshold'))

        elif self.algorithm == 'IMAT':
            if not all(elem in kwargs.keys() for elem in ['core', 'tolerance', 'epsilon']):
                warnings.warn('IMAT properties missing. The following arguments must be present: core, '
                              'tolerance, epsilon')
                return None
            else:
                # Get indexes of core reactions' ids:
                core_idx = [self.generic_model.reactions.index(rxn) for rxn in kwargs.get('core')]
                # Create IMATProperties object:
                return IMATProperties(exp_vector=list(self.reaction_calls.values()), exp_thresholds=[1, 2],
                                      core=core_idx, tolerance=kwargs.get('tolerance'), epsilon=kwargs.get('epsilon'))

        elif self.algorithm == 'tINIT':
            if not all(elem in kwargs.keys() for elem in
                       ['present_metabolites', 'essential_reactions', 'production_weight', 'allow_excretion',
                        'no_reverse_loops', 'solver']):
                warnings.warn('tINIT properties missing. The following arguments must be present: present_metabolites, '
                              'essential_reactions, production_weight, allow_excretion, no_reverse_loops, solver')
                return None
            else:
                # Get indexes of present_metabolites' ids:
                present_metabolites = kwargs.get('present_metabolites')
                if present_metabolites is not None:
                    present_metabolites_idx = [self.generic_model.metabolites.index(rxn)
                                               for rxn in present_metabolites]
                else:
                    present_metabolites_idx = None
                # Get indexes of essential_ractions' ids:
                essential_reactions = kwargs.get('essential_reactions')
                if essential_reactions is not None:
                    essential_reactions_idx = [self.generic_model.reactions.index(rxn)
                                               for rxn in essential_reactions]
                else:
                    essential_reactions_idx = None
                # Create tINITProperties object:
                return tINITProperties(reactions_scores=list(self.reaction_calls.values()),
                                       present_metabolites=present_metabolites_idx,
                                       essential_reactions=essential_reactions_idx,
                                       production_weight=kwargs.get('production_weight'),
                                       allow_excretion=kwargs.get('allow_excretion'),
                                       no_reverse_loops=kwargs.get('no_reverse_loops'), solver=kwargs.get('solver'))

        elif self.algorithm == 'CORDA':
            if not all(elem in kwargs.keys() for elem in ['pr_to_np', 'constraint', 'constrainby', 'om', 'ntimes', 'nl',
                                                          'solver']):
                warnings.warn('CORDA properties missing. The following arguments must be present: pr_to_np, '
                              'constraint, constrainby, om, ntimes, nl, solver')
                return None
            else:
                # Get High confidence reactions (ids) and indexes of said reactions:
                high_rxns_id = [k for (k, v) in self.reaction_calls.items() if v == 'High']
                high_rxns_idx = [self.generic_model.reactions.index(rxn) for rxn in high_rxns_id]
                # Get Medium confidence reactions (ids) and indexes oreconstructed_modelf said reactions:
                medium_rxns_id = [k for (k, v) in self.reaction_calls.items() if v == 'Medium']
                medium_rxns_idx = [self.generic_model.reactions.index(rxn) for rxn in medium_rxns_id]
                # Get Negative confidence reactions (ids) and indexes of said reactions:
                negative_rxns_id = [k for (k, v) in self.reaction_calls.items() if v == 'Negative']
                negative_rxns_idx = [self.generic_model.reactions.index(rxn) for rxn in negative_rxns_id]
                # Create CORDAProperties object:
                return CORDAProperties(high_rxns_idx, medium_rxns_idx, negative_rxns_idx,
                                       pr_to_np=kwargs.get('pr_to_np'), constraint=kwargs.get('pr_to_np'),
                                       constrainby=kwargs.get('constrainby'), om=kwargs.get('om'),
                                       ntimes=kwargs.get('ntimes'), nl=kwargs.get('nl'), solver=kwargs.get('solver'))

    def run(self):
        x = ModelBasedWrapper(self.generic_model)
        y = ReconstructionWrapper(x)
        self.reconstructed_model = y.run(self.algorithm_properties)
