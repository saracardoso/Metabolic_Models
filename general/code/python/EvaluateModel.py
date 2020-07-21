from os.path import join
from os import getcwd
from pandas import read_csv, DataFrame
import numpy as np

from troppo.tasks.task_io import JSONTaskIO
from troppo.tasks.core import TaskEvaluator


class EvaluateModel(object):

    def __init__(self, model):
        self.model = model

        # Get the tasks to evaluate:
        tasks_file = join(getcwd(), 'general/utility_data/metabolic_tasks_all_cells.json')
        task_reader = JSONTaskIO()
        self.tasks = task_reader.read_task(tasks_file)
        self.tasks_result = None

        # Get the different media to consider:
        self.media = read_csv(join(getcwd(), 'general/utility_data/media.csv'), index_col='ID').to_dict()
        self.media_result = None

    def evaluate_tasks(self):
        res = np.zeros((len(self.tasks), 2))
        with self.model as model_to_test:
            for boundary in model_to_test.boundary:
                boundary.knock_out()
            task_evaluator = TaskEvaluator(model=model_to_test, tasks=self.tasks)
            for i, task in enumerate(task_evaluator.tasks):
                task_evaluator.current_task = task
                res[i, ] = [task, task_evaluator.evaluate()]
        self.tasks_result = res[:]

    def evaluate_media_biomass_capacity(self):
        reaction_ids = [i.id for i in self.model.reactions]
        biomass_reactions = [i for i in reaction_ids if i.startswith('BIOMASS_')]

        res = DataFrame(np.zeros((len(self.media)+1, len(biomass_reactions))),
                        columns=biomass_reactions, index=self.media.keys())
        with self.model as model_to_test:
            for i_m, medium in enumerate(self.media.keys()):
                model_to_test.medium = self.media[medium]
                for i_r, biomass_reaction in enumerate(biomass_reactions):
                    model_to_test.objective = biomass_reaction
                    res.loc[i_m, i_r] = model_to_test.slim_optimize()
        self.media_result = res[:]
