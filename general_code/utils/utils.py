from cobra.io import read_sbml_model, write_sbml_model
from cobra.flux_analysis.variability import find_blocked_reactions
from os.path import join
from os import getcwd

# Base Directories:
base_dir = getcwd()

models_dir = join(base_dir, '0Models')
recon3d_model_dir = join(models_dir, 'Recon3D')

general_code_dir = join(base_dir, 'general_code')
general_code_utils_dir = join(general_code_dir, 'utils')
utils_entrez_genes_dir = join(general_code_utils_dir, 'entrez_genes')
utils_GPRs_dir = join(general_code_utils_dir, 'GPRs')

'''
Recon3D
'''
# Read recon3d original model:
recon3d_original = read_sbml_model(join(recon3d_model_dir, 'Recon3D.xml.gz'))
# Get model's genes:
gene_numbers = []
for i in recon3d_original.genes: gene_numbers.append(i.id)
f = open(join(utils_entrez_genes_dir, 'recon3D_genes.txt'), 'w')
for item in gene_numbers: f.write("%s\n" % item)
f.close()
# Get model's GPRs:
reaction_ids = []
for i in recon3d_original.reactions: reaction_ids.append(i.id)
f = open(join(utils_GPRs_dir, 'recon3D_GPR.txt'), 'w')
for item in reaction_ids: f.write("%s\t%s\n" % (item, recon3d_original.reactions.get_by_id(item).gene_reaction_rule))
f.close()



'''
Recon3D_consistent
'''
recon3d_consistent = recon3d_original.copy()
# Remove sink (sink_) and drug module (DM_) reactions:
recon3d_consistent.remove_reactions(
    [r for r in recon3d_consistent.reactions if r.id[:3] == 'DM_' or r.id[:5] == 'sink_'], remove_orphans=True)
# Remove reactions that became blocked after the removal of these reactions:
blocked_reactions = find_blocked_reactions(recon3d_consistent)
recon3d_consistent.remove_reactions(blocked_reactions, remove_orphans=True)
# Save consistent version of Recon3D into a file:
write_sbml_model(recon3d_consistent, join(recon3d_model_dir, 'Recon3D_consistent.xml.gz'))
# Get model's genes:
gene_numbers = []
for i in recon3d_consistent.genes: gene_numbers.append(i.id)
f = open(join(utils_entrez_genes_dir, 'recon3D_consistent_genes.txt'), 'w')
for item in gene_numbers: f.write("%s\n" % item)
f.close()
# Get model's GPRs:
reaction_ids = []
for i in recon3d_consistent.reactions: reaction_ids.append(i.id)
f = open(join(utils_GPRs_dir, 'recon3D_consistent_GPR.txt'), 'w')
for item in reaction_ids: f.write("%s\t%s\n" % (item, recon3d_consistent.reactions.get_by_id(item).gene_reaction_rule))
f.close()



'''
Recon3D_forTcells
'''
from cobra import Reaction
recon3d_for_Tcells = recon3d_consistent.copy()
# Create a reaction for the biomass reaction from the macrophage model iAB-AM0-1410:
new_reaction = Reaction(id='biomass_mac', name='Biomass reaction from macrophage model iAB-AM0-1410',
                        subsystem='Exchange/demand reaction',
                        lower_bound=0., upper_bound=1000.)
recon3d_for_Tcells.add_reactions([new_reaction])
recon3d_for_Tcells.reactions.biomass_mac.add_metabolites({
    'ala_L[c]': -0.396559456,
    'alpa_hs[c]': -0.011499127,
    'amp[c]': -0.048664064,
    'arg_L[c]': -0.325724532,
    'asn_L[c]': -0.215407845,
    'asp_L[c]': -0.282759085,
    'atp[c]': -25.17352552,
    'chsterol[c]': -0.020930954,
    'cmp[c]': -0.042373167,
    'cys_L[c]': -0.127154496,
    'dag_hs[c]': -0.0036682,
    'damp[c]': -0.021495345,
    'dcmp[c]': -0.014937443,
    'dgmp[c]': -0.014937443,
    'dtmp[c]': -0.021495345,
    'gln_L[c]': -0.280436629,
    'glu_L[c]': -0.424428935,
    'gly[c]': -0.366948135,
    'glygn1[c]': -0.528027894,
    'gmp[c]': -0.043710887,
    'h2o[c]': -25.17352552,
    'hdca[c]': -0.004850777,
    'hdcea[c]': -0.001222285,
    'his_L[c]': -0.153862747,
    'ile_L[c]': -0.25953452,
    'leu_L[c]': -0.580614138,
    'lys_L[c]': -0.351852168,
    'met_L[c]': -0.126573882,
    'ocdca[c]': -0.004736708,
    'ocdcea[c]': -0.003853116,
    'pail_hs[c]': -0.003741686,
    'pchol_hs[c]': -0.031527146,
    'pe_hs[c]': -0.021107135,
    'pglyc_hs[c]': -0.008918017,
    'phe_L[c]': -0.214246617,
    'pro_L[c]': -0.346626641,
    'ps_hs[c]': -0.001024655,
    'ser_L[c]': -0.476684207,
    'sphmyln_hs[c]': -0.007049706,
    'tag_hs[c]': -0.002742439,
    'thr_L[c]': -0.303661194,
    'trp_L[c]': -0.069673697,
    'ttdca[c]': -0.00136164,
    'tyr_L[c]': -0.156185203,
    'ump[c]': -0.04602478,
    'val_L[c]': -0.347207255,
    'adp[c]': 25.17352552,
    'h[c]': 25.17352552,
    'pi[c]': 25.17352552
})
recon3d_for_Tcells.reactions.biomass_mac.gene_reaction_rule = ''
# Save this version of Recon3D into a file:
write_sbml_model(recon3d_consistent, join(recon3d_model_dir, 'Recon3D_forTcells.xml.gz'))
# Get model's genes:
gene_numbers = []
for i in recon3d_consistent.genes: gene_numbers.append(i.id)
f = open(join(utils_entrez_genes_dir, 'recon3D_forTcells_genes.txt'), 'w')
for item in gene_numbers: f.write("%s\n" % item)
f.close()
# Get model's GPRs:
reaction_ids = []
for i in recon3d_consistent.reactions: reaction_ids.append(i.id)
f = open(join(utils_GPRs_dir, 'recon3D_forTcells_GPR.txt'), 'w')
for item in reaction_ids: f.write("%s\t%s\n" % (item, recon3d_consistent.reactions.get_by_id(item).gene_reaction_rule))
f.close()