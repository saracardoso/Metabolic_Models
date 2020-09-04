from pandas import read_csv, DataFrame
import numpy as np
import warnings

from troppo.tasks.task_io import JSONTaskIO
from troppo.tasks.core import TaskEvaluator


class EvaluateModel(object):

    def __init__(self, model, tasks_file=None, media_file=None):
        self.model = model.copy()

        # Get the tasks to evaluate:
        if tasks_file is not None:
            task_reader = JSONTaskIO()
            self.tasks = task_reader.read_task(tasks_file)
        else:
            self.tasks = None
        self.tasks_result = None

        # Get the different media to consider:
        if media_file is not None:
            self.get_media_from_file(media_file)
        else:
            self.media = None
        self.media_result = None

    def get_tasks_from_file(self, tasks_file):
        task_reader = JSONTaskIO()
        self.tasks = task_reader.read_task(tasks_file)
        self.tasks_result = None

    def get_media_from_file(self, media_file):
        self.media = read_csv(media_file, index_col='ID').to_dict()
        self.media['Default'] = ''
        self.media_result = None

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

        reaction_ids = [i.id for i in self.model.reactions]
        biomass_reactions = [i for i in reaction_ids if i.startswith('biomass_')]

        res = DataFrame(np.zeros((len(self.media), len(biomass_reactions))),
                        columns=biomass_reactions, index=self.media.keys())
        with self.model as model_to_test:
            default_medium = model_to_test.medium.copy()
            for i_m, medium in enumerate(self.media.keys()):
                if medium == 'Default':
                    model_to_test.medium = default_medium.copy()
                else:
                    model_to_test.medium = self.media[medium].copy()
                for i_r, biomass_reaction in enumerate(biomass_reactions):
                    model_to_test.objective = biomass_reaction
                    res.iloc[i_m, i_r] = model_to_test.slim_optimize()
        self.media_result = res[:]

    def save_tasks_result_csv(self, file_name):
        self.tasks_result[0].to_csv(file_name)

    def save_media_biomass_capacity_csv(self, file_name):
        self.media_result.to_csv(file_name)
