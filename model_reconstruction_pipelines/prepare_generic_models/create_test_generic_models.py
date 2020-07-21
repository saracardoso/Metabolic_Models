from os.path import join
from os import getcwd

from cobra.io import read_sbml_model, write_sbml_model
from cobra.manipulation.delete import remove_genes
from cobra.flux_analysis.variability import find_blocked_reactions
from cobra import Reaction

from general_code.python.EvaluateModel import EvaluateModel

# Base Directories:
base_dir = getcwd()
models_dir = join(base_dir, '0Models')
recon3d_model_dir = join(models_dir, 'Recon3D')

# Read recon3d original model:
recon3d = read_sbml_model(join(recon3d_model_dir, 'Recon3D.xml.gz'))

'''
Recon3d_forCancer
'''


# ---Get Model---

recon3d_forCancer = recon3d.copy()
# Remove transcript 8041.1 (and 0?)
remove_genes(recon3d_forCancer, ['8041_AT1', '0'], remove_reactions=False)
# Remove 'BIOMASS_maintenance' and 'BIOMASS_maintenance_noTrTr'
recon3d_forCancer.remove_reactions(['BIOMASS_maintenance', 'BIOMASS_maintenance_noTrTr'], remove_orphans=True)
# Remove artificial sinks and drug modules
recon3d_forCancer.remove_reactions(
    [r for r in recon3d_forCancer.reactions if r.id[:3] == 'DM_' or r.id[:3] == 'SK_'], remove_orphans=True)
# Remove resulting block reactions and respective genes and metabolites that end up with no reaction associated
blocked_reactions = find_blocked_reactions(recon3d_forCancer)
recon3d_forCancer.remove_reactions(blocked_reactions, remove_orphans=True)

# ---Evaluate Model---

evaluate_recon3d_forCancer = EvaluateModel(model=recon3d_forCancer)
# Task evaluation:
evaluate_recon3d_forCancer.evaluate_tasks()
# Test capacity to carry flux over 'biomass_reaction' using different mediums
# (recon3D medium,  Plasmax, Plasmax_v2, Plasmax_unconstrained, ):
# Save into data.frame the medium used and the flux through the biomass reaction
evaluate_recon3d_forCancer.evaluate_media_biomass_capacity()

# ---Write the model into a sbml file---
write_sbml_model(recon3d_forCancer, join(recon3d_model_dir, 'recon3d_forCancer.xml.gz'))


'''
Recon3d_forTcells
'''

#---Get Model---

recon3d_forTcells = recon3d.copy()
# Remove transcript 8041.1 (and 0?)
recon3d_forTcells = remove_genes(recon3d_forTcells, ['8041_AT1', '0'], remove_reactions=False)
# Remove 'BIOMASS_maintenance_noTrTr'
recon3d_forTcells.remove_reactions(['BIOMASS_maintenance_noTrTr'], remove_orphans=True)
# Add the biomass reaction from the macrophage model iAB-AM0-1410:
new_reaction = Reaction(id='biomass_mac', name='Biomass reaction from macrophage model iAB-AM0-1410',
                        subsystem='Exchange/demand reaction',
                        lower_bound=0., upper_bound=1000.)
recon3d_forTcells.add_reactions([new_reaction])
recon3d_forTcells.reactions.biomass_mac.add_metabolites({
    'ala_L[c]': -0.396559456, 'alpa_hs[c]': -0.011499127, 'amp[c]': -0.048664064, 'arg_L[c]': -0.325724532,
    'asn_L[c]': -0.215407845, 'asp_L[c]': -0.282759085, 'atp[c]': -25.17352552,
    'chsterol[c]': -0.020930954, 'cmp[c]': -0.042373167, 'cys_L[c]': -0.127154496,
    'dag_hs[c]': -0.0036682, 'damp[c]': -0.021495345, 'dcmp[c]': -0.014937443, 'dgmp[c]': -0.014937443,
    'dtmp[c]': -0.021495345,
    'gln_L[c]': -0.280436629, 'glu_L[c]': -0.424428935, 'gly[c]': -0.366948135, 'glygn1[c]': -0.528027894,
    'gmp[c]': -0.043710887,
    'h2o[c]': -25.17352552, 'hdca[c]': -0.004850777, 'hdcea[c]': -0.001222285, 'his_L[c]': -0.153862747,
    'ile_L[c]': -0.25953452,
    'leu_L[c]': -0.580614138, 'lys_L[c]': -0.351852168,
    'met_L[c]': -0.126573882,
    'ocdca[c]': -0.004736708, 'ocdcea[c]': -0.003853116,
    'pail_hs[c]': -0.003741686, 'pchol_hs[c]': -0.031527146, 'pe_hs[c]': -0.021107135, 'pglyc_hs[c]': -0.008918017,
    'phe_L[c]': -0.214246617, 'pro_L[c]': -0.346626641, 'ps_hs[c]': -0.001024655,
    'ser_L[c]': -0.476684207, 'sphmyln_hs[c]': -0.007049706,
    'tag_hs[c]': -0.002742439, 'thr_L[c]': -0.303661194, 'trp_L[c]': -0.069673697, 'ttdca[c]': -0.00136164,
    'tyr_L[c]': -0.156185203,
    'ump[c]': -0.04602478,
    'val_L[c]': -0.347207255,
    'adp[c]': 25.17352552,
    'h[c]': 25.17352552,
    'pi[c]': 25.17352552
})
recon3d_forTcells.reactions.biomass_mac.gene_reaction_rule = ''
# Remove artificial sinks and drug modules
recon3d_forTcells.remove_reactions(
    [r for r in recon3d_forTcells.reactions if r.id[:3] == 'DM_' or r.id[:5] == 'sink_'], remove_orphans=True)
# Remove resulting block reactions and respective genes and metabolites that end up with no reaction associated
blocked_reactions = find_blocked_reactions(recon3d_forTcells)
recon3d_forTcells.remove_reactions(blocked_reactions, remove_orphans=True)
# Write the model io a sbml file:
write_sbml_model(recon3d_forTcells, join(recon3d_model_dir, 'recon3d_forTcells.xml.gz'))

#---Evaluate Model---








