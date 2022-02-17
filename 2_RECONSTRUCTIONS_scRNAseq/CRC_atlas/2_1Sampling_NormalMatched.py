if __name__ == '__main__':
    from os.path import join, isdir
    from os import listdir, mkdir
    from dill import load
    import json
    from numpy import nan
    from gc import collect
    from pandas import read_csv

    from GENERAL.code.python.Validation import sampling

    '''
    Directories and files
    '''

    # Directories:
    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'
    utilityData_dir = join(base_dir, 'GENERAL/utility_data')
    CRCReconstructionNormalMatched_dir = join(base_dir, '2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched')

    # Files:
    tcells_media_file = join(utilityData_dir, 'media_Tcells_8percSerum.csv')
    sampling_evalTcells_file = join(utilityData_dir, 'sampling_evals_Tcells.json')

    '''
    Tcells
    '''

    # 0. Tcells models that will be sampled:
    abTcells = ['Naive CD4 Tcells', 'Memory CD4 Tcells', 'Proliferative CD4 Tcells', 'Regulatory CD4 Tcells',
                'IL17+ CD4 Tcells', 'Follicular CD4 Tcells', 'Naive CD8 Tcells', 'Cytotoxic CD8 Tcells',
                'Memory CD8 Tcells', 'Proliferative CD8 Tcells']

    # 1. Normal media, no additional constraints:
    individuals = listdir(CRCReconstructionNormalMatched_dir)
    for indiv in individuals:
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
                    print('--', cell_type)
                    if cell_type in ['Naive CD4 Tcells', 'Proliferative CD4 Tcells', 'Regulatory CD4 Tcells',
                                     'IL17+ CD4 Tcells', 'Follicular CD4 Tcells', 'Naive CD8 Tcells',
                                     'Memory CD8 Tcells', 'Proliferative CD8 Tcells']: # abTcells:
                        biomass_res = read_csv(join(CRCReconstructionNormalMatched_dir, indiv,
                                                    ''.join(('3_biomassAfter_', samp_name, '.csv'))),
                                               header=0, index_col=0)
                        if biomass_res.loc[cell_type, 'Blood_SMDB'] == 0:
                            print('This cell-type model will not be sampled: no biomass under Blood_SMDB medium.')
                        else:
                            sampling_result = sampling(model, tcells_media_file, change_media_bounds=None,
                                                       change_internal_bounds=None, n_samples=1000, processes=5)
                            ct_dir = cell_type.replace(' ', '_')
                            sampling_dir = join(CRCReconstructionNormalMatched_dir, indiv, '1_sampling')
                            if not isdir(sampling_dir):
                                mkdir(sampling_dir)
                                mkdir(join(sampling_dir, 'control'))
                            if not isdir(join(sampling_dir, 'control', samp_name)):
                                mkdir(join(sampling_dir, 'control', samp_name))
                            mkdir(join(sampling_dir, 'control', samp_name, ct_dir))
                            for medium, medium_result in sampling_result.items():
                                if medium_result is not None:
                                    medium_result.to_csv(join(sampling_dir, 'control', samp_name, ct_dir,
                                                              ''.join((medium, '.csv'))))
                    else:
                        print('This cell-type model will not be sampled')
                del temp_dump
                collect()
'''
    # 2. Changed media and/or additional constraints:
    with open(sampling_evalTcells_file, 'r') as f:
        sampling_evals_Tcells = json.load(f)
    individuals = listdir(CRCReconstructionNormalMatched_dir)
    for indiv in individuals:
        print('\n', indiv)
        # Get files that starts with 02_
        indiv_samples = [file for file in listdir(join(CRCReconstructionNormalMatched_dir, indiv))
                         if file.startswith('02_')]
        for samp in indiv_samples:
            samp_name = samp.replace('.csv', '').replace('02_', '')
            print('- ', samp_name)
            with open(join(CRCReconstructionNormalMatched_dir, indiv, samp), 'rb') as dump_file:
                temp_dump = load(dump_file)
                for cell_type, model in temp_dump.items():
                    print('--', cell_type)
                    if cell_type in abTcells:
                        ct_dir = join(sampling_dir, 'control', samp_name, cell_type.replace(' ', '_'))
                        sampling_files = listdir(ct_dir)
                        samp_results = {}
                        for samp_file in sampling_files:
                            samp_results[samp_file.replace('.csv', '')] = read_csv(join(ct_dir, samp_file),
                                                                                   header=0, index_col=0)
                        for evalID, sampling_evals in sampling_evals_Tcells.items():
                            print('---', evalID)
                            new_change_internal_bounds = None
                            if sampling_evals['change_internal_bounds'] is not None:
                                new_change_internal_bounds = {}
                                for medium, medium_results in samp_results.items():
                                    new_change_internal_bounds[medium] = {}
                                    for k, i in sampling_evals['change_internal_bounds'].items():
                                        if i[0] is nan and i[1] is nan:
                                            m_val = min(abs(medium_results.loc[:, k]))
                                            i[0] = -m_val
                                            i[1] = m_val
                                        elif i[0] is nan:
                                            i[0] = max(medium_results.loc[:, k])
                                        else:
                                            sampling_values = medium_results.loc[:, k]
                                            i[1] = min(sampling_values[sampling_values > 0])
                                    new_change_internal_bounds[medium][k] = i
                            sampling_result = sampling(model, tcells_media_file,
                                                       change_media_bounds=sampling_evals['change_media_bounds'],
                                                       change_internal_bounds=new_change_internal_bounds,
                                                       n_samples=1000, processes=6)
                            if not isdir(join(CRCReconstructionNormalMatched_dir, indiv, '1_sampling', 'not_control')):
                                mkdir(join(sampling_dir, 'not_control'))
                                mkdir(join(sampling_dir, 'not_control', samp_name))
                            mkdir(join(sampling_dir, 'not_control', samp_name, ct_dir))
                            for medium, medium_result in sampling_result.items():
                                medium_result.to_csv(join(sampling_dir, 'not_control', samp_name, ct_dir,
                                                          ''.join((medium, '.csv'))))
                    else:
                        print('This cell-type model will not be sampled')
                del temp_dump
                collect()

########################################################################################################################

# TCELLS

# GET MEDIA FOR TCELLS
from pandas import read_csv

media = read_csv('/home/scardoso/Documents/PhD/Metabolic_Models/GENERAL/utility_data/media_8percSerum.csv',
                 index_col='ID')

# SAMPLING FOR IN VIVO VS IN VITRO

# MAS PRIMEIRO, FAZER UM EXEMPLO DE SAMPLING COM OPEN BOUNDS NAS NAIVE cd4 OLHAR PARA A BIOMASS, OXPHOS AND GLUTAMINOLYSIS
from dill import load

with open(
        '/home/scardoso/Documents/PhD/Metabolic_Models/2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/31/02_scrEXT001.obj',
        'rb') as file:
    models_31 = load(file)

from cobra.sampling import OptGPSampler

import matplotlib.pyplot as plt
from datetime import datetime

print('Initializing sampler...')
sampler = OptGPSampler(models_31.reconstructed_models['Naive CD4 Tcells'], processes=4)
print('Sampling...')
s1 = sampler.sample(100)
print('Done!')

s1_plot = plt.figure()
plt.hist(s1.loc[s1.loc[:, 'MAR13082'] != 0, 'MAR13082'], range=(-0.5, 20), bins=205,
         histtype='stepfilled', edgecolor='#002800', facecolor='#E2F0CB')
plt.scatter(0, sum(s1.loc[:, 'MAR13082'] == 0), c='red', marker='X', linewidths=3)
plt.text(0, sum(s1.loc[:, 'MAR13082'] == 0), '(0,{})'.format(sum(s1.loc[:, 'MAR13082'] == 0)), color='red')
plt.xlabel('Flux (mmol gDW-1 hr-1)')
plt.ylabel('Nº Samples')
plt.title('Biomass (MAR13082)')
s1_plot.show()

print('Initializing sampler...')
print('Start: ', datetime.now().strftime('%H:%M:%S'))
sampler2 = OptGPSampler(models_31.reconstructed_models['Proliferative CD4 Tcells'], processes=4)
print('Sampling...')
s2 = sampler2.sample(100)
print('End: ', datetime.now().strftime('%H:%M:%S'))
print('Done!')

s2_plot = plt.figure()
plt.hist(s2.loc[s2.loc[:, 'MAR13082'] != 0, 'MAR13082'], range=(-0.5, 20), bins=205,
         histtype='stepfilled', edgecolor='#002800', facecolor='#E2F0CB')
plt.scatter(0, sum(s2.loc[:, 'MAR13082'] == 0), c='red', marker='X', linewidths=3)
plt.text(0, sum(s2.loc[:, 'MAR13082'] == 0), '(0,{})'.format(sum(s2.loc[:, 'MAR13082'] == 0)), color='red')
plt.xlabel('Flux (mmol gDW-1 hr-1)')
plt.ylabel('Nº Samples')
plt.title('Biomass (MAR13082)')
s2_plot.show()

with models_31.reconstructed_models['Naive CD4 Tcells'] as model:
    model.medium = media['Plasmax_serum'].to_dict()
    print('Initializing sampler...')
    print('Start: ', datetime.now().strftime('%H:%M:%S'))
    sampler3 = OptGPSampler(model, processes=4)
    print('Sampling...')
    s3 = sampler3.sample(100)
    print('End: ', datetime.now().strftime('%H:%M:%S'))
    print('Done!')

s3_plot = plt.figure()
plt.hist(s3.loc[s3.loc[:, 'MAR13082'] != 0, 'MAR13082'], range=(-0.5, 20), bins=205,
         histtype='stepfilled', edgecolor='#002800', facecolor='#E2F0CB')
plt.scatter(0, sum(s3.loc[:, 'MAR13082'] == 0), c='red', marker='X', linewidths=3)
plt.text(0, sum(s3.loc[:, 'MAR13082'] == 0), '(0,{})'.format(sum(s3.loc[:, 'MAR13082'] == 0)), color='red')
plt.xlabel('Flux (mmol gDW-1 hr-1)')
plt.ylabel('Nº Samples')
plt.title('Biomass (MAR13082)')
s3_plot.show()

# Check, using sampling, for expected phenotypes (organize into functions? and put them into the evaluation object?)
'''
