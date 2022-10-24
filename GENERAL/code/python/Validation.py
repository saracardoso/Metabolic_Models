from datetime import datetime
from os.path import join
from os import listdir

import numpy as np
from dill import load
from pandas import read_csv, concat, Series
from gc import collect
from re import search
from math import isnan

from cobra.sampling import OptGPSampler
from cobra.flux_analysis.room import room
from cobra.flux_analysis import pfba
from cobra import Solution


'''
Sampling
'''


def sampling(model, media_file, change_media_bounds: dict = None, change_internal_bounds: dict = None,
             n_samples: int = 1000, processes: int = 4):
    # if number in a change_media_bounds' value is a string, it will be multiplied by the value in the medium
    media = read_csv(media_file, index_col='ID')
    sampling_results = {}
    with model as model_test:
        for medium in media.columns:
            print('\n', datetime.now().strftime('%H:%M:%S'), ' | ', 'MEDIUM: ', medium)
            model_test.medium = media[medium].to_dict()
            if change_media_bounds is not None:
                print(datetime.now().strftime('%H:%M:%S'), ' | ', '- Changing media bounds..')
                new_media = model_test.medium.copy()
                for cpd_in, val in change_media_bounds[medium].items():
                    if isinstance(val, str): new_media[cpd_in] = float(val) * new_media[cpd_in]
                    else: new_media[cpd_in] = val
                model_test.medium = new_media.copy()
            if change_internal_bounds is not None:
                print(datetime.now().strftime('%H:%M:%S'), ' | ', '- Changing internal bounds..')
                for rxn, rxn_bounds in change_internal_bounds.values():
                    model_test.reactions.get_by_id(rxn).bounds = rxn_bounds
            print(datetime.now().strftime('%H:%M:%S'), ' | ', 'Initializing sampler...')
            sampler = OptGPSampler(model, processes=processes, thinning=100)
            # sampler = ACHRSampler(model, thinning=100)
            print(datetime.now().strftime('%H:%M:%S'), ' | ', 'Sampling...')
            try:
                fluxes_samples = sampler.sample(n_samples)
            except RuntimeError:
                sampling_results[medium] = None
                print('Did not work (Can not escape sampling region, model seems numerically unstable).')
            else:
                print(datetime.now().strftime('%H:%M:%S'), ' | ', 'Getting validation status...')
                val_status = sampler.validate(fluxes_samples)
                print(datetime.now().strftime('%H:%M:%S'), ' | ', 'Done!')
                print(datetime.now().strftime('%H:%M:%S'), ' | ', 'Filtering invalid samples...')
                sampling_results[medium] = fluxes_samples[val_status == 'v']
            print(datetime.now().strftime('%H:%M:%S'), ' | ', 'Done!')
    return sampling_results


'''
FBA
'''


def fba_evaluations(model, fva_results, fba_result, change_media_bounds: dict = None, change_internal_bounds: dict = None):
    # if number in a change_media_bounds' value is a string, it will be multiplied by the value in the medium
    with model as model_test:
        if change_media_bounds is not None:
            print('- Changing media bounds..')
            new_media = model_test.medium.copy()
            for cpd_in, val in change_media_bounds.items():
                if isinstance(val, str): new_media[cpd_in] = float(val) * new_media[cpd_in]
                else: new_media[cpd_in] = val
            model_test.medium = new_media.copy()
        count = 0
        if change_internal_bounds is not None:
            print('- Changing internal bounds..')
            new_change_internal_bounds = {}
            for rxn, rxn_bounds in change_internal_bounds.items():
                if isnan(rxn_bounds[0]) and isnan(rxn_bounds[1]):
                    m_val = min(abs(fva_results.loc[rxn, :]))
                    rxn_bounds[0] = -m_val
                    rxn_bounds[1] = m_val
                elif isnan(rxn_bounds[0]):
                    rxn_bounds[0] = fva_results.loc[rxn, 'maximum']
                elif isnan(rxn_bounds[1]):
                    rxn_bounds[1] = fva_results.loc[rxn, 'minimum']
                else:
                    count += 1
                    pass
                new_change_internal_bounds[rxn] = rxn_bounds
            for rxn, rxn_bounds in new_change_internal_bounds.items():
                model_test.reactions.get_by_id(rxn).bounds = rxn_bounds
        print(datetime.now().strftime('-Starting FBA...'))
        if count == 0:
            try:
                #solFBA = model_test.optimize(objective_sense='maximize')
                solFBA = pfba(model_test)
            except:
                print('Infeasible')
                solFBA = Series(np.zeros(11894).tolist())
            else:
                solFBA = solFBA.fluxes
        else:
            # ROOM
            solFBA = room(model_test, solution=fba_result, linear=True)
            solFBA = solFBA.fluxes
        print('-Done!')
    return solFBA#.fluxes


def helper_fba_condition_normalmatched(condition, media_file, medium, models_use, models_dir,
                                       individuals, fba_results=None, fva_results=None):
    predicted_fluxes = None
    models_names = []
    for indiv in individuals:
        print('\n', indiv)
        # Get files that starts with 02_
        indiv_samples = [file for file in listdir(join(models_dir, indiv)) if file.startswith('02_')]
        for samp in indiv_samples:
            samp_name = samp.replace('.obj', '').replace('02_', '')
            print('- ', samp_name)
            with open(join(models_dir, indiv, samp), 'rb') as dump_file:
                temp_dump = load(dump_file)
                for cell_type, model in temp_dump.items():
                    model_name = '_'.join((indiv, samp_name, cell_type))
                    if model_name in models_use:
                        print('--', cell_type)
                        # Get SMDB medium
                        media = read_csv(media_file, index_col='ID')
                        model.medium = media[medium].to_dict()
                        # Get objective
                        # ---
                        #fba_orignal_fluxes = fba_results.loc[:, model_name]
                        #if cell_type in ['Proliferative CD4 Tcells', 'Proliferative CD8 Tcells']:
                        #    model.objective = {model.reactions.MAR13082: 1}
                        #    fba_orig_value = fba_orignal_fluxes['MAR13082']
                        #else:
                        #    model.objective = {model.reactions.MAR13082: 1, model.reactions.MAR06916: 1}
                        #    fba_orig_value = fba_orignal_fluxes['MAR13082'] + fba_orignal_fluxes['MAR06916']
                        #fba_result = Solution(fba_orig_value, 'optimal', fluxes=fba_orignal_fluxes)
                        model.objective = {model.reactions.MAR13082: 1}
                        # ---
                        # Run FBA
                        #fluxes = fba_evaluations(model, fva_results['_'.join((indiv, samp_name, cell_type))],
                        #                         fba_results.loc[:,'_'.join((indiv, samp_name, cell_type))],
                        #                         change_media_bounds=condition['change_media_bounds'],
                        #                         change_internal_bounds=condition['change_internal_bounds'])
                        fluxes = fba_evaluations(model, None, None, #fba_result,
                                                 change_media_bounds=condition['change_media_bounds'],
                                                 change_internal_bounds=condition['change_internal_bounds'])
                        # Save predicted fluxes
                        models_names.append(model_name)
                        if predicted_fluxes is None:
                            predicted_fluxes = concat([fluxes], axis=1)
                        else:
                            predicted_fluxes = concat([predicted_fluxes, fluxes], axis=1)
                del temp_dump
                collect()
    predicted_fluxes.columns = models_names
    return predicted_fluxes
