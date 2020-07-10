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
Recon3D_Model_301
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
Recon3D_Model_301_consistent
'''
# Create consistent model and save into a file:
recon3d_consistent = recon3d_original.copy()
recon3d_consistent.remove_reactions(
    [r for r in recon3d_consistent.reactions if r.id[:3] == 'DM_' or r.id[:5] == 'sink_'], remove_orphans=True)
blocked_reactions = find_blocked_reactions(recon3d_consistent)
recon3d_consistent.remove_reactions(blocked_reactions, remove_orphans=True)
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
