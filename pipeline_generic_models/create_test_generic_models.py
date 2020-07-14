from os.path import join
from os import getcwd

from cobra.io import read_sbml_model, write_sbml_model
from cobra.flux_analysis.variability import find_blocked_reactions

# Base Directories:
base_dir = getcwd()

models_dir = join(base_dir, '0Models')
recon3d_model_dir = join(models_dir, 'Recon3D')

# Read recon3d original model:
recon3d = read_sbml_model(join(recon3d_model_dir, 'Recon3D.xml.gz'))

'''
Recon3d_forCancer
'''


#---Get Model---

recon3d_forCancer = recon3d.copy()
# Remove transcript 8041.1 (and 0?)

# Remove 'biomass_maintenaince' and 'biomass_maintenance_noTrTr'

# Remove artificial sinks and drug modules

# Remove resulting block reactions and respective genes and metabolites that end up with no reaction associated


#---Evaluate Model---

#












