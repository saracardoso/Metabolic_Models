# -*- coding: utf-8 -*-

import logging
import multiprocessing
from builtins import dict, map
from functools import partial
from itertools import product
from typing import List, Set, Union

import pandas as pd
from optlang.exceptions import SolverError

from cobra.core import Configuration, Gene, Reaction
from cobra.flux_analysis.moma import add_moma
from cobra.flux_analysis.room import add_room
from cobra.manipulation.delete import find_gene_knockout_reactions
from cobra.util import solver as sutil


LOGGER = logging.getLogger(__name__)


CONFIGURATION = Configuration()


def _reactions_knockouts_with_restore(model, reactions):
    with model:
        for reaction in reactions:
            reaction.knock_out()
        growth = _get_growth(model)
    return [r.id for r in reactions], growth, model.solver.status


def _get_growth(model):
    try:
        if "moma_old_objective" in model.solver.variables:
            model.slim_optimize()
            growth = model.solver.variables.moma_old_objective.primal
        else:
            growth = model.slim_optimize()
    except SolverError:
        growth = float("nan")
    return growth


def _reaction_deletion(model, ids):
    return _reactions_knockouts_with_restore(
        model, [model.reactions.get_by_id(r_id) for r_id in ids]

    )


def _gene_deletion(model, ids):
    all_reactions = []
    for g_id in ids:
        all_reactions.extend(find_gene_knockout_reactions(model, (model.genes.get_by_id(g_id),)))
    _, growth, status = _reactions_knockouts_with_restore(model, all_reactions)
    return (ids, growth, status)


def _reaction_deletion_worker(ids):
    global _model
    return _reaction_deletion(_model, ids)


def _gene_deletion_worker(ids):
    global _model
    return _gene_deletion(_model, ids)


def _init_worker(model):
    global _model
    _model = model


def _multi_deletion(model, entity, element_lists, method="fba", solution=None, processes=None, **kwargs):
    """
    Provide a common interface for single or multiple knockouts.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to perform deletions in.
    entity : 'gene' or 'reaction'
        The entity to knockout (``cobra.Gene`` or ``cobra.Reaction``).
    element_lists : list
        List of iterables ``cobra.Reaction``s or ``cobra.Gene``s (or their IDs)
        to be deleted.
    method: {"fba", "moma", "linear moma", "room", "linear room"}, optional
        Method used to predict the growth rate.
    solution : cobra.Solution, optional
        A previous solution to use as a reference for (linear) MOMA or ROOM.
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to the number of CPUs found.
    kwargs :
        Passed on to underlying simulation functions.

    Returns
    -------
    pandas.DataFrame
        A representation of all combinations of entity deletions. The
        columns are 'growth' and 'status', where

        index : tuple(str)
            The gene or reaction identifiers that were knocked out.
        growth : float
            The growth rate of the adjusted model.
        status : str
            The solution's status.
    """
    solver = sutil.interface_to_str(model.problem.__name__)
    if method == "moma" and solver not in sutil.qp_solvers:
        raise RuntimeError(
            "Cannot use MOMA since '{}' is not QP-capable."
            "Please choose a different solver or use FBA only.".format(solver)
        )

    if processes is None:
        processes = CONFIGURATION.processes

    with model:
        if "moma" in method:
            add_moma(model, solution=solution, linear="linear" in method)
        elif "room" in method:
            add_room(model, solution=solution, linear="linear" in method, **kwargs)

        args = set([frozenset(comb) for comb in product(*element_lists)])
        processes = min(processes, len(args))

        def extract_knockout_results(result_iter):
            result = pd.DataFrame(
                [
                    (
                        set(ids),
                        growth,
                        status,
                    )
                    for (ids, growth, status) in result_iter
                ],
                columns=["ids", "growth", "status"],
            )
            return result

        if processes > 1:
            worker = dict(gene=_gene_deletion_worker, reaction=_reaction_deletion_worker)[entity]
            chunk_size = len(args) // processes
            pool = multiprocessing.Pool(processes, initializer=_init_worker, initargs=(model,))
            results = extract_knockout_results(pool.imap_unordered(worker, args, chunksize=chunk_size))
            pool.close()
            pool.join()
        else:
            worker = dict(gene=_gene_deletion, reaction=_reaction_deletion)[entity]
            results = extract_knockout_results(map(partial(worker, model), args))
        return results


def _entities_ids(entities):
    try:
        return [e.id for e in entities]
    except AttributeError:
        return list(entities)


def _element_lists(entities, *ids):
    lists = list(ids)
    if lists[0] is None:
        lists[0] = entities
    result = [_entities_ids(lists[0])]
    for l in lists[1:]:
        if l is None:
            result.append(result[-1])
        else:
            result.append(_entities_ids(l))
    return result


def single_gene_deletion(model, gene_list=None, method="fba", solution=None, processes=None, **kwargs):
    """
    Knock out each gene from a given list.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to perform deletions in.
    gene_list : iterable
        ``cobra.Gene``s to be deleted. If not passed,
        all the genes from the model are used.
    method: {"fba", "moma", "linear moma", "room", "linear room"}, optional
        Method used to predict the growth rate.
    solution : cobra.Solution, optional
        A previous solution to use as a reference for (linear) MOMA or ROOM.
    processes : int, optional
        The number of parallel processes to run. Can speed up the computations
        if the number of knockouts to perform is large. If not passed,
        will be set to the number of CPUs found.
    kwargs :
        Keyword arguments are passed on to underlying simulation functions
        such as ``add_room``.

    Returns
    -------
    pandas.DataFrame
        A representation of all single gene deletions. The columns are
        'growth' and 'status', where

        index : tuple(str)
            The gene identifier that was knocked out.
        growth : float
            The growth rate of the adjusted model.
        status : str
            The solution's status.

    """
    return _multi_deletion(
        model,
        "gene",
        element_lists=_element_lists(model.genes, gene_list),
        method=method,
        solution=solution,
        processes=processes,
        **kwargs)


#######################################


if __name__ == '__main__':
    from os.path import join
    from os import listdir
    from dill import load
    import json
    from pandas import read_csv, concat, DataFrame
    from gc import collect
    from re import search
    from numpy import arange

    from cobra import Solution
    from cobra.flux_analysis.room import room

    '''
    Directories and files
    '''

    # Directories:
    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'
    utilityData_dir = join(base_dir, 'GENERAL/utility_data')
    CRCReconstruction_dir = join(base_dir, '2_RECONSTRUCTIONS_scRNAseq/CRC_atlas')
    CRCReconstructionNormalMatched_dir = join(CRCReconstruction_dir, 'NormalMatched')
    HumanGEM_dir = join(base_dir, '0MODELS/HumanGEM')

    # Files:
    tcells_media_file = join(utilityData_dir, 'media_Tcells_8percSerum.csv')
    single_gene_deletion_file = join(utilityData_dir, 'single_gene_deletion.json')
    CRCatlas_sampling_file = join(CRCReconstruction_dir, 'CRCatlas_sampling.json')

    '''
    Run single gene deletions
    '''

    CRCatlas_sampling = json.load(open(CRCatlas_sampling_file))
    genes_rxns_toTest = json.load(open(single_gene_deletion_file))
    predicted_fluxes = read_csv(join(CRCReconstructionNormalMatched_dir, '1_control_analysis/FBA/normal_FBA.csv'),
                                index_col=0)
    individuals = listdir(CRCReconstructionNormalMatched_dir)
    predicted_fluxes2 = None
    objective_values = DataFrame(0, index=list(genes_rxns_toTest.keys()), columns=arange(172))
    models_names = []
    for gene_id, rxns in genes_rxns_toTest.items():
        print('\n', gene_id)
        models_names = []
        predicted_fluxes2 = None
        obj_vals = []
        for indiv in individuals:
            if search('[.]', indiv):
                next
            print('\n', indiv)
            # Get files that starts with 02_
            indiv_samples = [file for file in listdir(join(CRCReconstructionNormalMatched_dir, indiv))
                             if file.startswith('02_')]
            for samp in indiv_samples:
                samp_name = samp.replace('.obj', '').replace('02_', '')
                print('- ', samp_name)
                with open(join(CRCReconstructionNormalMatched_dir, indiv, samp), 'rb') as dump_file:
                    temp_dump = load(dump_file)
                    for cell_type, model in temp_dump.items():
                        if cell_type in CRCatlas_sampling['NormalMatched'][indiv]['control'][samp_name]:
                            print('--', cell_type)
                            # Get SMDB medium
                            media = read_csv(tcells_media_file, index_col='ID')
                            model.medium = media['Blood_SMDB'].to_dict()
                            # Get objective
                            fba_orignal_fluxes = predicted_fluxes.loc[:, '_'.join((indiv, samp_name, cell_type))]
                            if cell_type in ['Proliferative CD4 Tcells', 'Proliferative CD8 Tcells']:
                                model.objective = {model.reactions.MAR13082: 1}
                                fba_orig_value = fba_orignal_fluxes['MAR13082']
                            else:
                                model.objective = {model.reactions.MAR13082: 1, model.reactions.MAR06916: 1}
                                fba_orig_value = fba_orignal_fluxes['MAR13082'] + fba_orignal_fluxes['MAR06916']
                            # Get original FBA without knockouts:
                            fba_result = Solution(fba_orig_value, 'optimal', fluxes=fba_orignal_fluxes)
                            # Get new predictions.
                            for rxn in rxns:
                                model.reactions.get_by_id(rxn).bounds = (0, 0)
                            sol = room(model, solution=fba_result, linear=True)
                            # Save predicted objective value:
                            obj_vals.append(sol.objective_value)
                            # Save predicted fluxes
                            models_names.append('_'.join((indiv, samp_name, cell_type)))
                            if predicted_fluxes2 is None:
                                predicted_fluxes2 = concat([sol.fluxes], axis=1)
                            else:
                                predicted_fluxes2 = concat([predicted_fluxes2, sol.fluxes], axis=1)
                del temp_dump
                collect()
        objective_values.loc[gene_id, :] = obj_vals
        predicted_fluxes2.columns = models_names
        predicted_fluxes2.to_csv(''.join((CRCReconstructionNormalMatched_dir,
                                          '/3_gene_essentiality/specific_genes/', gene_id, '.csv')))
    objective_values.columns = models_names
    objective_values.to_csv(''.join((CRCReconstructionNormalMatched_dir,
                                     '/3_gene_essentiality/objective_values.csv')))

    '''
    Check glucose case:
    '''

    to_test = {}
    to_test['rxn'] = {}
    to_test['rxn']['SLC5A1'] = ['MAR01378', 'MAR06895', 'MAR08884']
    to_test['rxn']['GLUT1'] = ['MAR05029']
    to_test['rxn']['MAR05450'] = ['MAR05450']
    to_test['medium'] = {}
    to_test['medium']['glucose'] = ['MAR09034']
    to_test['medium']['glucose_sucrose'] = ['MAR09034', 'MAR09416']
    to_test['medium']['glucose_gcpool'] = ['MAR09034', 'MAR12043']
    to_test['medium']['glucose_sucrose_gcpool'] = ['MAR09034', 'MAR12043']

    CRCatlas_sampling = json.load(open(CRCatlas_sampling_file))
    predicted_fluxes = read_csv(join(CRCReconstructionNormalMatched_dir, '1_control_analysis/FBA/normal_FBA.csv'),
                                index_col=0)
    individuals = listdir(CRCReconstructionNormalMatched_dir)
    for test_type in to_test.keys():
        for test_id in to_test[test_type].keys():
            print('\n', test_type)
            models_names = []
            predicted_fluxes2 = None
            for indiv in individuals:
                if search('[.]', indiv):
                    next
                print('\n', indiv)
                # Get files that starts with 02_
                indiv_samples = [file for file in listdir(join(CRCReconstructionNormalMatched_dir, indiv))
                                 if file.startswith('02_')]
                for samp in indiv_samples:
                    samp_name = samp.replace('.obj', '').replace('02_', '')
                    print('- ', samp_name)
                    with open(join(CRCReconstructionNormalMatched_dir, indiv, samp), 'rb') as dump_file:
                        temp_dump = load(dump_file)
                        for cell_type, model in temp_dump.items():
                            if cell_type in CRCatlas_sampling['NormalMatched'][indiv]['control'][samp_name]:
                                print('--', cell_type)
                                # Get SMDB medium
                                media = read_csv(tcells_media_file, index_col='ID')
                                model.medium = media['Blood_SMDB'].to_dict()
                                # Get objective
                                fba_orignal_fluxes = predicted_fluxes.loc[:, '_'.join((indiv, samp_name, cell_type))]
                                if cell_type in ['Proliferative CD4 Tcells', 'Proliferative CD8 Tcells']:
                                    model.objective = {model.reactions.MAR13082: 1}
                                    fba_orig_value = fba_orignal_fluxes['MAR13082']
                                else:
                                    model.objective = {model.reactions.MAR13082: 1, model.reactions.MAR06916: 1}
                                    fba_orig_value = fba_orignal_fluxes['MAR13082'] + fba_orignal_fluxes['MAR06916']
                                # Get original FBA without knockouts:
                                fba_result = Solution(fba_orig_value, 'optimal', fluxes=fba_orignal_fluxes)
                                # Get new predictions.
                                if test_type == 'rxn':
                                    with model as model_test:
                                        for rxn in to_test[test_type][test_id]:
                                            model_test.reactions.get_by_id(rxn).bounds = (0, 0)
                                        sol = room(model_test, solution=fba_result, linear=True)
                                else:
                                    with model as model_test:
                                        new_media = model_test.medium.copy()
                                        for metab in to_test[test_type][test_id]:
                                            new_media[metab] = 0
                                        model_test.medium = new_media.copy()
                                        sol = model_test.optimize(objective_sense='maximize')
                                # Save predicted fluxes
                                models_names.append('_'.join((indiv, samp_name, cell_type)))
                                if predicted_fluxes2 is None:
                                    predicted_fluxes2 = concat([sol.fluxes], axis=1)
                                else:
                                    predicted_fluxes2 = concat([predicted_fluxes2, sol.fluxes], axis=1)
                    del temp_dump
                    collect()
            predicted_fluxes2.columns = models_names
            predicted_fluxes2.to_csv(''.join((CRCReconstructionNormalMatched_dir,
                                              '/3_gene_essentiality/glucose_case/', test_id, '.csv')))
