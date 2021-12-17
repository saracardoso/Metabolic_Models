from pandas import read_csv, DataFrame
import numpy as np
import warnings

from troppo.tasks.task_io import JSONTaskIO
from troppo.tasks.core import TaskEvaluator

from troppo.methods_wrappers import GapfillWrapper
from troppo.methods.gapfill.efm import DEFAULT_CONFIG

from cobra.flux_analysis import find_blocked_reactions


class EvaluateModel(object):

    def __init__(self, model, tasks_file=None, media_file=None,
                 biomass_reactions=None, biomass_metab='MAM03971s'):
        self.model = model.copy()

        # Get the tasks to evaluate:
        if tasks_file is not None:
            self._get_tasks_from_file(tasks_file)
        else:
            self.tasks = None
        self.tasks_result = None

        # Get the different media to consider:
        if media_file is not None:
            self._get_media_from_file(media_file)
        else:
            self.media = None
        self.media_result = None
        if biomass_reactions is None:
            self.biomass_reactions = ['MAR13082']
        else:
            self.biomass_reactions = biomass_reactions
        self.biomass_metabolite = biomass_metab

    def _get_tasks_from_file(self, tasks_file):
        task_reader = JSONTaskIO()
        self.tasks = task_reader.read_task(tasks_file)

    def _get_media_from_file(self, media_file):
        self.media = read_csv(media_file, index_col='ID').to_dict()
        self.media['Default'] = ''

    def evaluate_tasks(self):
        if self.tasks is None:
            warnings.warn('No tasks available.')
            return

        tasks_ids = [i.name for i in self.tasks]
        res_df = DataFrame(np.zeros((len(self.tasks), 4)),
                           columns=['Task_Name', 'Task_Group', 'Observed', 'Expected'],
                           index=tasks_ids)
        res_dict = dict()
        with self.model as model_to_test:
            for boundary in model_to_test.boundary:
                boundary.knock_out()
            task_evaluator = TaskEvaluator(model=model_to_test, tasks=self.tasks, solver='CPLEX')
            for i, task in enumerate(task_evaluator.tasks):
                task_evaluator.current_task = task
                task_id = self.tasks[i].name
                task_name = self.tasks[i].annotations['name']
                task_group = self.tasks[i].annotations['task_group']
                expected = not self.tasks[i].should_fail

                sol = task_evaluator.evaluate()
                # sol gives True if should fail and failed, and if should not fail and did not fail
                # sol will be changed so that everytime it works is True and everytime it fails is False:
                if not expected:
                    observed = not sol[0]
                else:
                    observed = sol[0]

                print(task_id, '|\t', task_name, '|\t', task_group, '|\t', observed, '-', expected)
                res_df.loc[task_id] = [task_name] + [task_group] + [observed] + [expected]
                res_dict[task_id] = sol[1].var_values()
        self.tasks_result = (res_df, res_dict)

    def evaluate_media_biomass_capacity(self):
        if self.media is None:
            warnings.warn('No media dataframe.')
            return

        res = DataFrame(np.zeros((len(self.media), len(self.biomass_reactions))),
                        columns=self.biomass_reactions, index=self.media.keys())
        with self.model as model_to_test:
            default_medium = model_to_test.medium.copy()
            for i_m, medium in enumerate(self.media.keys()):
                if medium == 'Default':
                    model_to_test.medium = {key: val for key, val in default_medium.copy().items() if val != 0}
                else:
                    model_to_test.medium = self.media[medium].copy()
                for i_r, biomass_reaction in enumerate(self.biomass_reactions):
                    model_to_test.objective = biomass_reaction
                    res.iloc[i_m, i_r] = model_to_test.slim_optimize()
        self.media_result = res[:]

    def save_tasks_result_csv(self, file_name):
        self.tasks_result[0].to_csv(file_name)

    def save_media_biomass_capacity_csv(self, file_name):
        self.media_result.to_csv(file_name)


class GapFillModel(EvaluateModel):

    def __init__(self, model, generic_model, tasks_file=None, media_file=None):
        super().__init__(model, tasks_file, media_file)
        self.generic_model = generic_model.copy()
        self.history = {0: model.copy()}
        self.sol_gapFill_tasks = dict()
        self.sol_gapFill_media = dict()

    def _run_gapfill_media(self):
        with self.model as model_to_test:
            default_medium = model_to_test.medium.copy()
            # Reactions not in model:
            not_present = []
            for rxn in model_to_test.reactions:
                if rxn.bounds == (0, 0):
                    not_present = not_present + [rxn.id]
            # Gap-fill media, if they failed:
            for medium in self.media_result.index:
                for biomass in self.media_result.columns:
                    if self.media_result.loc[medium, biomass] == 0:
                        print('Gap filling medium ', medium, ' for biomass reaction ', biomass)
                        if medium != 'Default':
                            model_to_test.medium = {key: val for key, val in self.media[medium].copy().items()
                                                    if val != 0}
                        else:
                            model_to_test.medium = default_medium
                        model_to_test.objective = biomass
                        gp = GapfillWrapper(model=model_to_test)
                        DEFAULT_CONFIG['BIG_M_VALUE'] = 1e6
                        sol = gp.run(avbl_fluxes=not_present,
                                     ls_override={'produced': [self.biomass_metabolite]},
                                     algorithm='efm',
                                     kshproperties=DEFAULT_CONFIG)
                        if len(sol) > 0:
                            self.sol_gapFill_media[medium + '_' + biomass] = sol[0]
                        else:
                            self.sol_gapFill_media[medium + '_' + biomass] = []

    def _add_reactions(self):
        history_new_key = max(self.history.keys()) + 1
        self.history[history_new_key] = self.model.copy()
        reactions_add = []
        for rxns in self.sol_gapFill_media.values():
            reactions_add = reactions_add + rxns
        for rxn in np.unique(reactions_add):
            self.model.reactions.get_by_id(rxn).bounds = self.generic_model.reactions.get_by_id(rxn).bounds

    def run(self):
        if self.media_result is None:
            print('\nThe capacity to produce biomass from the media given has not been yet assessed. '
                  'Will do that now...')
            self.evaluate_media_biomass_capacity()

        # Inspect if gapfill is at all needed:
        count = 0
        for medium in self.media_result.index:
            for biomass in self.media_result.columns:
                if self.media_result.loc[medium, biomass] == 0:
                    count += 1
        if count == 0:
            print('\nNo GapFill needed, all biomass reactions to test produce flux in all media provided! ')
        else:
            print('')
            self._run_gapfill_media()
            self._add_reactions()

    def final_model(self):
        print('\nInfo before removing blocked reactions:')
        self._print_model_info(self.model)
        new_model = self.model.copy()
        blocked_reactions = find_blocked_reactions(new_model, open_exchanges=True)
        new_model.remove_reactions(blocked_reactions, remove_orphans=True)
        print('\nInfo after removing blocked reactions (final model):')
        self._print_model_info(new_model)
        return new_model

    def _print_model_info(self, model):
        model_reactions = model.reactions
        model_genes = model.genes
        model_metabolites = model.metabolites
        print('- Number of Reactions:', len(model_reactions))
        print('  ... Number of those with bounds not (0,0):', len([rxn for rxn in model_reactions
                                                                   if rxn.bounds != (0, 0)]))
        print('- Number of Genes:', len(model_genes))
        print('- Number of Metabolites:', len(model_metabolites))
        pass
