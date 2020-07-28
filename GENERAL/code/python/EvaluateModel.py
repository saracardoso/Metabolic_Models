from os.path import join
from os import getcwd
from pandas import read_csv, DataFrame
import numpy as np
import warnings

from troppo.tasks.task_io import JSONTaskIO
from troppo.tasks.core import TaskEvaluator


class EvaluateModel(object):

    def __init__(self, model, media_file=None):
        self.model = model.copy()

        # Get the tasks to evaluate:
        tasks_file = join(getcwd(), 'general/utility_data/metabolic_tasks_all_cells.json')
        task_reader = JSONTaskIO()
        self.tasks = task_reader.read_task(tasks_file)
        self.tasks_result = None

        # Get the different media to consider:
        if media_file is not None:
            self.get_media_from_file(media_file)
        else:
            self.media = None
        self.media_result = None

    def get_media_from_file(self, media_file):
        self.media = read_csv(media_file, index_col='ID').to_dict()
        self.media['Default'] = ''

    def evaluate_tasks(self):
        tasks_descriptions = [i.annotations['description'] for i in self.tasks]
        res_df = DataFrame(np.zeros((len(self.tasks), 2)),
                           columns=['Observed', 'Expected'],
                           index=tasks_descriptions)
        res_dict = dict()
        with self.model as model_to_test:
            for boundary in model_to_test.boundary:
                boundary.knock_out()
            task_evaluator = TaskEvaluator(model=model_to_test, tasks=self.tasks, solver='CPLEX')
            for i, task in enumerate(task_evaluator.tasks):
                task_evaluator.current_task = task
                task_description = self.tasks[i].annotations['description']
                expected = not self.tasks[i].should_fail

                sol = task_evaluator.evaluate()
                res_df.loc[task_description] = [sol[0]] + [expected]
                res_dict[task_description] = sol[1]
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
                    model_to_test.medium = default_medium
                else:
                    model_to_test.medium = self.media[medium]
                for i_r, biomass_reaction in enumerate(biomass_reactions):
                    model_to_test.objective = biomass_reaction
                    res.iloc[i_m, i_r] = model_to_test.slim_optimize()
        self.media_result = res[:]
