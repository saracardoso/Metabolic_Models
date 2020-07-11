#Directories:
base_dir = getwd() # Open through project file

Omics_data_Repos = 'home/scardoso/Documents/PhD/Omics_Data_Repos'
T_Cell_Data_Repo = paste(Omics_data_Repos, 'T_Cell_Data_Repo', sep='/')
T_Cell_Data_Repo_Process_Scripts = paste(T_Cell_Data_Repo, 'Process_Scripts', sep='/')
T_Cell_Data_Repo_Process_Scripts_utils = paste(T_Cell_Data_Repo_Process_Scripts, 'utils', sep='/')

Metabolic_Models_Repo_general_code = paste(base_dir, 'general_code', sep='/')
Metabolic_Models_Repo_general_utils = paste(Metabolic_Models_Repo_general_code, 'utils', sep='/')
Metabolic_Models_Repo_general_utils_entrez_genes = paste(Metabolic_Models_Repo_general_utils, 'entrez_genes', sep='/')
Metabolic_Models_Repo_general_utils_GPRs = paste(Metabolic_Models_Repo_general_utils, 'GPRs', sep='/')

#Files:
genes_lengths_GRCh38p13_file = paste(T_Cell_Data_Repo_Process_Scripts_utils, 'genes_lengths_GRCh38p13.txt', sep='/')

recon_genes_orig_file = paste(Metabolic_Models_Repo_general_utils, 'recon_genes_orig.tsv', sep='/')

recon3DModel_301_genes_file = paste(Metabolic_Models_Repo_general_utils_entrez_genes, 'recon3D_genes.txt', sep='/')
recon3DModel_301_gene_mapping_new_file = paste(Metabolic_Models_Repo_general_utils, 'recon3D_gene_mapping.csv', sep='/')

recon3DModel_301_consistent_genes_file = paste(Metabolic_Models_Repo_general_utils_entrez_genes, 'recon3D_consistent_genes.txt', sep='/')

recon3D_GPR_file = paste(Metabolic_Models_Repo_general_utils_GPRs, 'recon3D_GPR.txt', sep='/')
recon3D_consistent_GPR_file = paste(Metabolic_Models_Repo_general_utils_GPRs, 'recon3D_consistent_GPR.txt', sep='/')

########################################
###MAP RECON3D GENES TO ENSEMBL GENES###
########################################


#Read file with all ensembl genes from our annotated genome (GRCh38.p13):
ensembl_genes = read.table(genes_lengths_GRCh38p13_file, header = T)$gene
#Read file with information on the genes present in the recon (entrez):
recon_genes = read.csv(recon_genes_orig_file, sep='\t')
recon_genes$gene_number = as.character(recon_genes$gene_number)
for(x in 1:length(recon_genes$gene_number)) recon_genes$gene_number[x] = strsplit(recon_genes$gene_number[x], '[.]')[[1]][1]
#Read file with genes (entrez) in the recon3d model:
recon3dmodel_gene_ids = as.character(read.table(recon3DModel_301_genes_file)$V1)

#Get unique recon genes:
unique_recon_genes = recon_genes[!duplicated(recon_genes$gene_number),]
View(unique_recon_genes)

#Get unique recon 3d model genes:
recon3dmodels_gene_ids_unique = c()
for (x in strsplit(recon3dmodel_gene_ids, '[_]')) recon3dmodels_gene_ids_unique = unique(c(recon3dmodels_gene_ids_unique, x[1]))

#Get info only for the recon genes in the recon3d model:
recon3dmodel_genes_info = unique_recon_genes[match(recon3dmodels_gene_ids_unique, unique_recon_genes$gene_number),]
recon3dmodel_genes_info = recon3dmodel_genes_info[-which(is.na(recon3dmodel_genes_info$gene_number)),]

#Get recon3d genes positions in our ensembl_genes:
map_data_genes_to_recon3d_genes = match(recon3dmodel_genes_info$ensembl_gene,ensembl_genes)

#First gather the recon genes that have ensembl info:
recon3d_entrez_ensembl = recon3dmodel_genes_info$ensembl_gene[!is.na(map_data_genes_to_recon3d_genes)]
names(recon3d_entrez_ensembl) = recon3dmodel_genes_info$gene_number[!is.na(map_data_genes_to_recon3d_genes)]

#For the recon genes that do not have ensembl info, add them manually:
genes_not_in_data = recon3dmodel_genes_info$gene_number[is.na(map_data_genes_to_recon3d_genes)]
#View(recon3dmodel_genes_info[match(genes_not_in_data, recon3dmodel_genes_info$gene_number),])
recon3d_entrez_ensembl = c(recon3d_entrez_ensembl,
                           'ENSG00000103222', 'ENSG00000175711', 'ENSG00000204310', 'ENSG00000012779', 'ENSG00000213760',
                           'ENSG00000235863', 'ENSG00000136881', 'ENSG00000161267', 'ENSG00000152254', 'ENSG00000073734',
                           'ENSG00000167703', 'ENSG00000142319', 'ENSG00000143554', 'ENSG00000204386', 'ENSG00000159131',
                           'ENSG00000177628', 'ENSG00000104522', 'ENSG00000277161', 'ENSG00000103415', 'ENSG00000204228',
                           'ENSG00000113504', 'ENSG00000167494', 'ENSG00000265788', 'ENSG00000231852', 'ENSG00000149527',
                           'ENSG00000276293', 'ENSG00000069764', 'ENSG00000162040', 'ENSG00000196620', 'ENSG00000221988',
                           'ENSG00000147813', 'ENSG00000157881', 'ENSG00000242612', 'ENSG00000224586', 'ENSG00000104524',
                           'ENSG00000226278', 'ENSG00000278540', 'ENSG00000006530', 'ENSG00000196139', 'ENSG00000198848',
                           'ENSG00000105852', 'ENSG00000073969', 'ENSG00000184983', 'ENSG00000170906', 'ENSG00000187581') #ensembl
names(recon3d_entrez_ensembl)[1838:1882] = genes_not_in_data #entrez
#'ENSG00000187581'%in%ensembl_genes
#Data frame with full info on recon 3d model genes:
recon3DModel_genes = data.frame(entrez=names(recon3d_entrez_ensembl),
                                ensembl=recon3d_entrez_ensembl,
                                symbol=recon3dmodel_genes_info$symbol[match(names(recon3d_entrez_ensembl),
                                                                            recon3dmodel_genes_info$gene_number)])
View(recon3DModel_genes)
#NOTE: There are two ensembl genes that are duplicated (each one is mapped to two different entrez ids)

write.csv(recon3DModel_genes, recon3DModel_301_gene_mapping_new_file, row.names=F)





#########################
##Get unique entrez ids##
#########################
recon3dmodel_consistent_gene_ids = as.character(read.table(recon3DModel_301_consistent_genes_file)$V1)
recon3dmodel_consistent_gene_ids_unique = c()
for (x in strsplit(recon3dmodel_consistent_gene_ids, '[_]')) recon3dmodel_consistent_gene_ids_unique = unique(c(recon3dmodel_consistent_gene_ids_unique, x[1]))
write.table(recon3dmodel_consistent_gene_ids_unique, recon3DModel_301_consistent_genes_file, row.names=F, col.names=F)

recon3dmodel_gene_ids = as.character(read.table(recon3DModel_301_genes_file)$V1)
recon3dmodel_gene_ids_unique = c()
for (x in strsplit(recon3dmodel_gene_ids, '[_]')) recon3dmodel_gene_ids_unique = unique(c(recon3dmodel_gene_ids_unique, x[1]))
write.table(recon3dmodel_gene_ids_unique, recon3DModel_301_genes_file, row.names=F, col.names=F)





################################
##Get GPRs without transcripts##
################################

# Recon3D:
recon3D_GPRs = read.table(recon3D_GPR_file, sep='\t')
for(gpr in 1:length(recon3D_GPRs$V2))
  recon3D_GPRs$V2[gpr] = paste(strsplit(recon3D_GPRs$V2[gpr], '_AT[1-9]+')[[1]], collapse='')
write.table(recon3D_GPRs, recon3D_GPR_file, row.names=F, col.names=F)

# Recon3D_consistent:
recon3D_consistent_GPRs = read.table(recon3D_consistent_GPR_file, sep='\t')
for(gpr in 1:length(recon3D_consistent_GPRs$V2))
  recon3D_consistent_GPRs$V2[gpr] = paste(strsplit(recon3D_consistent_GPRs$V2[gpr], '_AT[1-9]+')[[1]], collapse='')
write.table(recon3D_consistent_GPRs, recon3D_consistent_GPR_file, row.names=F, col.names=F)





