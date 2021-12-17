
# TCELLS

# GET MEDIA FOR TCELLS
from pandas import read_csv

media = read_csv('/home/scardoso/Documents/PhD/Metabolic_Models/GENERAL/utility_data/media_8percSerum.csv', index_col='ID')

# SAMPLING FOR IN VIVO VS IN VITRO

# MAS PRIMEIRO, FAZER UM EXEMPLO DE SAMPLING COM OPEN BOUNDS NAS NAIVE cd4 OLHAR PARA A BIOMASS, OXPHOS AND GLUTAMINOLYSIS
from dill import load

with open('/home/scardoso/Documents/PhD/Metabolic_Models/2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/31/02_scrEXT001.obj', 'rb') as file:
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

# Check, using sampling, for expected phenotypes (organize into functions? and put them into the evaluation object?)
