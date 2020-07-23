from os.path import join
from os import getcwd

from cobra.io import read_sbml_model, write_sbml_model
from cobra.flux_analysis.variability import find_blocked_reactions
from cobra import Reaction

from general.code.python.EvaluateModel import EvaluateModel

# Base Directories:
base_dir = getcwd()
models_dir = join(base_dir, '0Models')
HumanGEM_dir = join(models_dir, 'HumanGEM')


'''
HumanGEM_1.4
'''

# --- Read HumanGEM original model ---

HumanGEM = read_sbml_model(join(HumanGEM_dir, 'HumanGEM_1.4.xml.gz'))

# --- Evaluate Model ---

evaluate_recon3d_forCancer = EvaluateModel(model=HumanGEM_forCancer)
# Task evaluation:
evaluate_recon3d_forCancer.evaluate_tasks()
# Test capacity to carry flux over 'biomass_reaction' using different mediums
# (recon3D medium,  Plasmax, Plasmax_v2, Plasmax_unconstrained, ):
# Save into data.frame the medium used and the flux through the biomass reaction
evaluate_recon3d_forCancer.evaluate_media_biomass_capacity()

'''
HumanGEM_forCancer
'''


# --- Get Model ---

HumanGEM_forCancer = HumanGEM.copy()
# Remove resulting block reactions and respective genes and metabolites that end up with no reaction associated
blocked_reactions = find_blocked_reactions(HumanGEM_forCancer)
HumanGEM_forCancer.remove_reactions(blocked_reactions, remove_orphans=True)

# --- Evaluate Model ---

evaluate_recon3d_forCancer = EvaluateModel(model=HumanGEM_forCancer)
# Task evaluation:
evaluate_recon3d_forCancer.evaluate_tasks()
# Test capacity to carry flux over 'biomass_reaction' using different mediums
# (recon3D medium,  Plasmax, Plasmax_v2, Plasmax_unconstrained, ):
# Save into data.frame the medium used and the flux through the biomass reaction
evaluate_recon3d_forCancer.evaluate_media_biomass_capacity()

# --- Write the model into a sbml file ---
write_sbml_model(HumanGEM_forCancer, join(HumanGEM_dir, 'HumanGEM_forCancer.xml.gz'))


'''
HumanGEM_forTcells
'''

#--- Get Model ---

HumanGEM_forTcells = HumanGEM.copy()
# Add the biomass maintenance reaction from Recon3D:
new_reaction_bm = Reaction(id='biomass_maintenance_Recon3D',
                           name='Biomass maintenance reaction from Recon3D',
                           subsystem='Exchange/demand reaction',
                           lower_bound=0., upper_bound=1000.)
HumanGEM_forTcells.add_reactions([new_reaction_bm])
HumanGEM_forTcells.reactions.biomass_maintenance_Recon3D.add_metabolites({
    'm01307c': -0.50563, 'm01365c': -0.35926, 'm01369c': -0.27942, 'm01370c': -0.35261, 'm01371c': -20.7045,
    'm01450c': -0.020401, 'm01589c': -0.011658, 'm01623c': -0.039036, 'm01628c': -0.046571, 'm01968c': -0.27519,
    'm01974c': -0.38587, 'm01975c': -0.326, 'm01986c': -0.53889, 'm02034c': -0.036117, 'm02040c': -20.6508,
    'm02125c': -0.12641, 'm02184c': -0.28608, 'm02360c': -0.54554, 'm02426c': -0.59211, 'm02471c': -0.15302,
    'm02684c': -0.15446, 'm02685c': -0.055374, 'm02715c': -0.002914, 'm02724c': -0.25947, 'm02750c': -0.023315,
    'm02770c': -0.41248, 'm02808c': -0.005829, 'm02896c': -0.39253, 'm02908c': -0.017486, 'm02993c': -0.31269,
    'm03089c': -0.013306, 'm03101c': -0.15967, 'm03130c': -0.053446, 'm03135c': -0.35261,
    'm01285c': 20.6508, 'm02039c': 20.6508, 'm02751c': 20.6508, 'temp001c': 1
})
HumanGEM_forTcells.reactions.biomass_mac.gene_reaction_rule = ''
# Add the biomass reaction from the macrophage model iAB-AM0-1410:
new_reaction_mac = Reaction(id='biomass_macrophage_iABAM01410',
                            name='Biomass reaction from macrophage model iAB-AM0-1410',
                            subsystem='Exchange/demand reaction',
                            lower_bound=0., upper_bound=1000.)
HumanGEM_forTcells.add_reactions([new_reaction_mac])
HumanGEM_forTcells.reactions.biomass_macrophage_iABAM01410.add_metabolites({
    'm01307c': -0.396559456, 'alpa_hs_c': -0.011499127, 'm01334c': -0.048664064, 'm01365c': -0.325724532,
    'm01369c': -0.215407845, 'm01370c': -0.282759085, 'm01371c': -25.17352552, 'm01450c': -0.020930954,
    'm01590c': -0.042373167, 'm01628c': -0.127154496, 'm00240c': -0.0036682, 'm01639c': -0.021495345,
    'm01644c': -0.014937443, 'm01686c': -0.014937443, 'm01752c': -0.021495345, 'm01975c[c]': -0.280436629,
    'm01974c': -0.424428935, 'm01986c': -0.366948135, 'm01990c': -0.528027894, 'm02016c': -0.043710887,
    'm02040c': -25.17352552, 'm02674c': -0.004850777, 'm02675c': -0.001222285, 'm02125c': -0.153862747,
    'm02184c': -0.25953452, 'm02360c': -0.580614138, 'm02426c': -0.351852168, 'm02471c': -0.126573882,
    'm02938c': -0.004736708, 'm02646l': -0.003853116, 'm02750c': -0.003741686, 'm02684c': -0.031527146,
    'm02685c': -0.021107135, 'm02715c': -0.008918017, 'm02724c': -0.214246617, 'm02770c': -0.346626641,
    'm02808c': -0.001024655, 'm02896c': -0.476684207, 'm02908c': -0.007049706, 'm02959c': -0.002742439,
    'm02993c': -0.303661194, 'm03089c': -0.069673697, 'm02494c': -0.00136164, 'm03101c': -0.156185203,
    'm03114c': -0.04602478, 'm03135c': -0.347207255,
    'm01285c': 25.17352552, 'm02039c': 25.17352552, 'm02751c': 25.17352552, 'temp001c': 1
})
HumanGEM_forTcells.reactions.biomass_mac.gene_reaction_rule = ''
# Remove resulting block reactions and respective genes and metabolites that end up with no reaction associated
blocked_reactions = find_blocked_reactions(HumanGEM_forTcells)
HumanGEM_forTcells.remove_reactions(blocked_reactions, remove_orphans=True)
# Write the model io a sbml file:
write_sbml_model(HumanGEM_forTcells, join(HumanGEM_dir, 'recon3d_forTcells.xml.gz'))

# --- Evaluate Model ---








