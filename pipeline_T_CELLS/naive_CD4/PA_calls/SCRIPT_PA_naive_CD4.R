#Directories:
base_dir = getwd() # When oppening with ptoject
Omics_data_Repos = 'home/scardoso/Documents/PhD/Omics_Data_Repos'

Metabolic_Models_Repo_general_code = paste(base_dir, 'general_code', sep='/')
Metabolic_Models_Repo_general_code_utils = paste(Metabolic_Models_Repo_general_code, 'utils', sep='/')
Metabolic_Models_Repo_general_code_R = paste(Metabolic_Models_Repo_general_code, 'R', sep='/')

Metabolic_Models_Repo_reconstruction_process_T_CELLS = paste(base_dir, 'reconstruction_process_T_CELLS', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD4 = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS, 'naive_CD4', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD4_PA_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD4, 'PA_calls', sep='/')

T_Cell_Data_Repo = paste(Omics_data_Repos, 'T_Cell_Data_Repo', sep='/')
T_Cell_Data_Repo_Data_Files = paste(T_Cell_Data_Repo, 'Data_Files', sep='/')
T_Cell_Data_Repo_Data_Files_naive_CD4 = paste(T_Cell_Data_Repo_Data_Files, 'naive_CD4', sep='/')
T_Cell_Data_Repo_Data_Files_naive_CD4_RNA_seq = paste(T_Cell_Data_Repo_Data_Files_naive_CD4, 'RNA_seq', sep='/')
T_Cell_Data_Repo_Data_Files_naive_CD4_Microarray = paste(T_Cell_Data_Repo_Data_Files_naive_CD4, 'Microarray', sep='/')
T_Cell_Data_Repo_Data_Files_naive_CD4_Proteomics = paste(T_Cell_Data_Repo_Data_Files_naive_CD4, 'Proteomics', sep='/')

#Files:
PA_calls_file = paste(Metabolic_Models_Repo_general_code_R, 'PA_calls.R', sep='/')
recon3D_consistent_genes_file = paste(Metabolic_Models_Repo_general_code_utils, 'entrez_genes/recon3D_consistent_genes.txt', sep='/')
recon3D_consistent_GPR_file = paste(Metabolic_Models_Repo_general_code_utils, 'GPRs/recon3D_consistent_GPR.txt', sep='/')
recon3DModel_gene_mapping_file = paste(Metabolic_Models_Repo_general_code_utils, 'recon3DModel_gene_mapping.csv', sep='/')

naive_CD4_microarray_data_file = paste(T_Cell_Data_Repo_Data_Files_naive_CD4_Microarray, 'naive_CD4_microarray_data.csv', sep='/')
naive_CD4_rnaseq_data_file = paste(T_Cell_Data_Repo_Data_Files_naive_CD4_RNA_seq, 'naive_CD4_rnaseq_data.csv', sep='/')
naive_CD4_proteomics_data_file = paste(T_Cell_Data_Repo_Data_Files_naive_CD4_Proteomics, 'naive_CD4_proteomics_data.csv', sep='/')
naive_CD4_microarray_metadata_file = paste(T_Cell_Data_Repo_Data_Files_naive_CD4_Microarray, 'naive_CD4_microarray_metadata.csv', sep='/')
naive_CD4_rnaseq_metadata_file = paste(T_Cell_Data_Repo_Data_Files_naive_CD4_RNA_seq, 'naive_CD4_rnaseq_metadata.csv', sep='/')

#Source:
source(PA_calls_file)


#--- Read expression data for PA calls ---

recon3D_consistent_genes = as.character(read.table(recon3D_consistent_genes_file)[,1])

naive_CD4_expression_data = read_omics_data(recon3D_consistent_genes, recon3DModel_gene_mapping_file,
                                            transcriptomics_data_files = list(RNAseq=naive_CD4_rnaseq_data_file, Microarray=naive_CD4_microarray_data_file),
                                            transcriptomics_metadata_files = list(RNAseq=naive_CD4_rnaseq_metadata_file, Microarray=naive_CD4_microarray_metadata_file),
                                            proteomics_data_files = list(ImmProt=naive_CD4_proteomics_data_file),
                                            proteomics_metadata_files = NULL)


#--- Reaction Calls ---














#-- Analysis of PA Calls --

#Heatmap with rnaseq and microarray samples and the final gene calls:
heatmap_gene_calls(naive_CD4_gene_PA_calls, plot_tilte='Naive CD4 | PA Gene Calls',
                   sample_colors_legend_location=c(65,0), cell_colours_legend_location=c(65,10))

#Cluster PA gene calls and see how they cluster (per type of data? and/or GSE? are there any clear outliers or clusters?)
naive_CD4_data_to_plot = list()
naive_CD4_data_to_plot$data = cbind(naive_CD4_gene_PA_calls$calls_per_sample$Transcriptomics$RNAseq,
                                    naive_CD4_gene_PA_calls$calls_per_sample$Transcriptomics$Microarray,
                                    naive_CD4_gene_PA_calls$calls_per_sample$Proteomics$ImmProt)
naive_CD4_data_to_plot$metadata = data.frame(GSE=c(naive_CD4_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                                   naive_CD4_expression_data$Transcriptomics$Microarray$metadata$GSE,
                                                   c(rep('Steady_State', 4), rep('Activated', 4))),
                                             data_type=c(rep('RNAseq', dim(naive_CD4_expression_data$Transcriptomics$RNAseq$data)[2]),
                                                         rep('Microarray', dim(naive_CD4_expression_data$Transcriptomics$Microarray$data)[2]),
                                                         rep('ImmProt', dim(naive_CD4_expression_data$Proteomics$ImmProt$data)[2])),
                                             row.names = c(rownames(naive_CD4_expression_data$Transcriptomics$RNAseq$metadata),
                                                           rownames(naive_CD4_expression_data$Transcriptomics$Microarray$metadata),
                                                           colnames(naive_CD4_gene_PA_calls$calls_per_sample$Proteomics$ImmProt)),
                                             stringsAsFactors = T)
cluster_analysis(naive_CD4_data_to_plot, color_leafs='GSE', label_leafs='data_type') #Color by GSE dataset and sample names by data type

#How many genes with score 0, 1 and 2:
sum(na.omit(naive_CD4_gene_PA_calls$final_gene_calls$Final==0)) # 1190
sum(na.omit(naive_CD4_gene_PA_calls$final_gene_calls$Final==1)) # 310
sum(na.omit(naive_CD4_gene_PA_calls$final_gene_calls$Final==2)) # 382

#How many reactions with score 0, 1 and 2:
sum(na.omit(naive_CD4_PA_reactions==0)) # 3080
sum(na.omit(naive_CD4_PA_reactions==1)) # 1437
sum(na.omit(naive_CD4_PA_reactions==2)) # 1421

#GO analysis:
# Get the genes from all the reactions present:
naive_CD4_reactions = names(na.omit(naive_CD4_PA_reactions))[na.omit(naive_CD4_PA_reactions)!=0]
naive_CD4_genes_from_reactions=recon3D_gpr_rules_filtered[match(naive_CD4_reactions, recon3D_gpr_rules_filtered[,1]),2]
naive_CD4_genes_from_reactions = unique(unlist(strsplit(unlist(strsplit(gsub('[)]','', gsub('[(]','', na.omit(naive_CD4_genes_from_reactions))), ' or ')), ' and ')))
naive_CD4_genes_from_reactions = naive_CD4_gene_PA_calls$final_gene_calls[match(naive_CD4_genes_from_reactions, rownames(naive_CD4_gene_PA_calls$final_gene_calls)),]
naive_CD4_genes_from_reactions = rownames(naive_CD4_genes_from_reactions)[naive_CD4_genes_from_reactions[,'Final']!=0]
naive_CD4_genes_from_reactions = recon3DModel_genes$ensembl[match(naive_CD4_genes_from_reactions, recon3DModel_genes$entrez)]
# Universe of genes: all genes in recon3D model:
universe_genes = unique(recon3DModel_genes$ensembl)
# GO analysis:
naive_CD4_GO_analysis_reactions = GO_analysis(naive_CD4_genes_from_reactions, universe_genes)
naive_CD4_GO_tables_reactions_sig = significant_GOs(naive_CD4_GO_analysis_reactions)
View(naive_CD4_GO_tables_reactions_sig$BP)
View(naive_CD4_GO_tables_reactions_sig$CC)
View(naive_CD4_GO_tables_reactions_sig$MF)







#For general GO:
# Get genes from reactions different between cell types and perform GO analysis.
# Hierarchical clustering using all samples, using final gene calls and using only the final reaction calls
# Do heatmap (colored with values from PA calls) using all samples, using final gene calls and using only final reaction calls

# How many genes in each GO are from each method, for those GO that are significant in at least one of the methods:
#m1_genes_in_go = c()
#m2_genes_in_go = c()
#gos = unique(c(rownames(M1_GO_tables_reactions_sig$BP),rownames(M2_GO_tables_reactions_sig$BP)))
#descriptions = unique(c(M1_GO_tables_reactions_sig$BP$Description,M2_GO_tables_reactions_sig$BP$Description))
#for (go in gos){
#  if(go%in%rownames(M1_GO_table_BP_reactions)) m1_genes_in_go = c(m1_genes_in_go, M1_GO_table_BP_reactions[go,'Count'])
#  else m1_genes_in_go = c(m1_genes_in_go, 0)
#  if(go%in%rownames(M2_GO_table_BP_reactions)) m2_genes_in_go = c(m2_genes_in_go, M2_GO_table_BP_reactions[go,'Count'])
#  else m2_genes_in_go = c(m2_genes_in_go, 0)
#}
#genes_in_go = data.frame(description=descriptions, M1=m1_genes_in_go, M2=m2_genes_in_go,
#                         row.names=gos)
#diff_genes_in_go = genes_in_go[genes_in_go$M1!=genes_in_go$M2,]
#View(diff_genes_in_go[order(abs(diff_genes_in_go$M1-diff_genes_in_go$M2), decreasing = T),])
#View(genes_in_go[genes_in_go$M1==genes_in_go$M2,])


all_types=list()
all_types$data = as.matrix(cbind(naive_CD4_transcriptomics_PA_calls$calls_per_sample$RNAseq,naive_CD4_transcriptomics_PA_calls$calls_per_sample$Microarray,
                                 Th1_transcriptomics_PA_calls$calls_per_sample$RNAseq,Th1_transcriptomics_PA_calls$calls_per_sample$Microarray,
                                 Th2_transcriptomics_PA_calls$calls_per_sample$RNAseq,Th2_transcriptomics_PA_calls$calls_per_sample$Microarray,
                                 Th17_transcriptomics_PA_calls$calls_per_sample$RNAseq,Th17_transcriptomics_PA_calls$calls_per_sample$Microarray,
                                 Treg_transcriptomics_PA_calls$calls_per_sample$RNAseq,Treg_transcriptomics_PA_calls$calls_per_sample$Microarray,
                                 Tfollicular_transcriptomics_PA_calls$calls_per_sample$RNAseq,Tfollicular_transcriptomics_PA_calls$calls_per_sample$Microarray,
                                 CD4_CM_transcriptomics_PA_calls$calls_per_sample$RNAseq,CD4_CM_transcriptomics_PA_calls$calls_per_sample$Microarray,
                                 CD4_EM_transcriptomics_PA_calls$calls_per_sample$RNAseq,CD4_EM_transcriptomics_PA_calls$calls_per_sample$Microarray,
                                 CD4_memory_transcriptomics_PA_calls$calls_per_sample$Microarray,
                                 naive_CD8_transcriptomics_PA_calls$calls_per_sample$RNAseq,naive_CD8_transcriptomics_PA_calls$calls_per_sample$Microarray,
                                 CD8_CM_transcriptomics_PA_calls$calls_per_sample$RNAseq,CD8_CM_transcriptomics_PA_calls$calls_per_sample$Microarray,
                                 CD8_EM_transcriptomics_PA_calls$calls_per_sample$RNAseq,CD8_EM_transcriptomics_PA_calls$calls_per_sample$Microarray
                                 ))
all_types$metadata = data.frame(GSE=c(naive_CD4_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                      naive_CD4_expression_data$Transcriptomics$Microarray$metadata$GSE,
                                      Th1_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                      Th1_expression_data$Transcriptomics$Microarray$metadata$GSE,
                                      Th2_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                      Th2_expression_data$Transcriptomics$Microarray$metadata$GSE,
                                      Th17_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                      Th17_expression_data$Transcriptomics$Microarray$metadata$GSE,
                                      Treg_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                      Treg_expression_data$Transcriptomics$Microarray$metadata$GSE,
                                      Tfollicular_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                      Tfollicular_expression_data$Transcriptomics$Microarray$metadata$GSE,
                                      CD4_CM_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                      CD4_CM_expression_data$Transcriptomics$Microarray$metadata$GSE,
                                      CD4_EM_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                      CD4_EM_expression_data$Transcriptomics$Microarray$metadata$GSE,
                                      CD4_memory_expression_data$Transcriptomics$Microarray$metadata$GSE,
                                      naive_CD8_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                      naive_CD8_expression_data$Transcriptomics$Microarray$metadata$GSE,
                                      CD8_CM_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                      CD8_CM_expression_data$Transcriptomics$Microarray$metadata$GSE,
                                      CD8_EM_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                      CD8_EM_expression_data$Transcriptomics$Microarray$metadata$GSE
                                      ),
                                data_type=c(rep('RNAseq', dim(naive_CD4_expression_data$Transcriptomics$RNAseq$data)[2]),
                                            rep('Microarray', dim(naive_CD4_expression_data$Transcriptomics$Microarray$data)[2]),
                                            rep('RNAseq', dim(Th1_expression_data$Transcriptomics$RNAseq$data)[2]),
                                            rep('Microarray', dim(Th1_expression_data$Transcriptomics$Microarray$data)[2]),
                                            rep('RNAseq', dim(Th2_expression_data$Transcriptomics$RNAseq$data)[2]),
                                            rep('Microarray', dim(Th2_expression_data$Transcriptomics$Microarray$data)[2]),
                                            rep('RNAseq', dim(Th17_expression_data$Transcriptomics$RNAseq$data)[2]),
                                            rep('Microarray', dim(Th17_expression_data$Transcriptomics$Microarray$data)[2]),
                                            rep('RNAseq', dim(Treg_expression_data$Transcriptomics$RNAseq$data)[2]),
                                            rep('Microarray', dim(Treg_expression_data$Transcriptomics$Microarray$data)[2]),
                                            rep('RNAseq', dim(Tfollicular_expression_data$Transcriptomics$RNAseq$data)[2]),
                                            rep('Microarray', dim(Tfollicular_expression_data$Transcriptomics$Microarray$data)[2]),
                                            rep('RNAseq', dim(CD4_CM_expression_data$Transcriptomics$RNAseq$data)[2]),
                                            rep('Microarray', dim(CD4_CM_expression_data$Transcriptomics$Microarray$data)[2]),
                                            rep('RNAseq', dim(CD4_EM_expression_data$Transcriptomics$RNAseq$data)[2]),
                                            rep('Microarray', dim(CD4_EM_expression_data$Transcriptomics$Microarray$data)[2]),
                                            rep('Microarray', dim(CD4_memory_expression_data$Transcriptomics$Microarray$data)[2]),
                                            rep('RNAseq', dim(naive_CD8_expression_data$Transcriptomics$RNAseq$data)[2]),
                                            rep('Microarray', dim(naive_CD8_expression_data$Transcriptomics$Microarray$data)[2]),
                                            rep('RNAseq', dim(CD8_CM_expression_data$Transcriptomics$RNAseq$data)[2]),
                                            rep('Microarray', dim(CD8_CM_expression_data$Transcriptomics$Microarray$data)[2]),
                                            rep('RNAseq', dim(CD8_EM_expression_data$Transcriptomics$RNAseq$data)[2]),
                                            rep('Microarray', dim(CD8_EM_expression_data$Transcriptomics$Microarray$data)[2])
                                            ),
                                cell_type = c(rep('naive_CD4', dim(naive_CD4_expression_data$Transcriptomics$RNAseq$data)[2]+dim(naive_CD4_expression_data$Transcriptomics$Microarray$data)[2]),
                                              rep('Th1', dim(Th1_expression_data$Transcriptomics$RNAseq$data)[2]+dim(Th1_expression_data$Transcriptomics$Microarray$data)[2]),
                                              rep('Th2', dim(Th2_expression_data$Transcriptomics$RNAseq$data)[2]+dim(Th2_expression_data$Transcriptomics$Microarray$data)[2]),
                                              rep('Th17', dim(Th17_expression_data$Transcriptomics$RNAseq$data)[2]+dim(Th17_expression_data$Transcriptomics$Microarray$data)[2]),
                                              rep('Treg', dim(Treg_expression_data$Transcriptomics$RNAseq$data)[2]+dim(Treg_expression_data$Transcriptomics$Microarray$data)[2]),
                                              rep('Tfollicular', dim(Tfollicular_expression_data$Transcriptomics$RNAseq$data)[2]+dim(Tfollicular_expression_data$Transcriptomics$Microarray$data)[2]),
                                              rep('CD4_CM', dim(CD4_CM_expression_data$Transcriptomics$RNAseq$data)[2]+dim(CD4_CM_expression_data$Transcriptomics$Microarray$data)[2]),
                                              rep('CD4_EM', dim(CD4_EM_expression_data$Transcriptomics$RNAseq$data)[2]+dim(CD4_EM_expression_data$Transcriptomics$Microarray$data)[2]),
                                              rep('CD4_memory', dim(CD4_memory_expression_data$Transcriptomics$Microarray$data)[2]),
                                              rep('naive_CD8', dim(naive_CD8_expression_data$Transcriptomics$RNAseq$data)[2]+dim(naive_CD8_expression_data$Transcriptomics$Microarray$data)[2]),
                                              rep('CD8_CM', dim(CD8_CM_expression_data$Transcriptomics$RNAseq$data)[2]+dim(CD8_CM_expression_data$Transcriptomics$Microarray$data)[2]),
                                              rep('CD8_EM', dim(CD8_EM_expression_data$Transcriptomics$RNAseq$data)[2]+dim(CD8_EM_expression_data$Transcriptomics$Microarray$data)[2])
                                              ),
                                row.names = c(rownames(naive_CD4_expression_data$Transcriptomics$RNAseq$metadata),
                                              rownames(naive_CD4_expression_data$Transcriptomics$Microarray$metadata),
                                              rownames(Th1_expression_data$Transcriptomics$RNAseq$metadata),
                                              rownames(Th1_expression_data$Transcriptomics$Microarray$metadata),
                                              rownames(Th2_expression_data$Transcriptomics$RNAseq$metadata),
                                              rownames(Th2_expression_data$Transcriptomics$Microarray$metadata),
                                              rownames(Th17_expression_data$Transcriptomics$RNAseq$metadata),
                                              rownames(Th17_expression_data$Transcriptomics$Microarray$metadata),
                                              rownames(Treg_expression_data$Transcriptomics$RNAseq$metadata),
                                              rownames(Treg_expression_data$Transcriptomics$Microarray$metadata),
                                              rownames(Tfollicular_expression_data$Transcriptomics$RNAseq$metadata),
                                              rownames(Tfollicular_expression_data$Transcriptomics$Microarray$metadata),
                                              rownames(CD4_CM_expression_data$Transcriptomics$RNAseq$metadata),
                                              rownames(CD4_CM_expression_data$Transcriptomics$Microarray$metadata),
                                              rownames(CD4_EM_expression_data$Transcriptomics$RNAseq$metadata),
                                              rownames(CD4_EM_expression_data$Transcriptomics$Microarray$metadata),
                                              rownames(CD4_memory_expression_data$Transcriptomics$RNAseq$metadata),
                                              rownames(CD4_memory_expression_data$Transcriptomics$Microarray$metadata),
                                              rownames(naive_CD8_expression_data$Transcriptomics$RNAseq$metadata),
                                              rownames(naive_CD8_expression_data$Transcriptomics$Microarray$metadata),
                                              rownames(CD8_CM_expression_data$Transcriptomics$RNAseq$metadata),
                                              rownames(CD8_CM_expression_data$Transcriptomics$Microarray$metadata),
                                              rownames(CD8_EM_expression_data$Transcriptomics$RNAseq$metadata),
                                              rownames(CD8_EM_expression_data$Transcriptomics$Microarray$metadata)
                                              ),
                                stringsAsFactors = T)
cluster_analysis(all_types, color_leafs='cell_type', label_leafs='GSE', leg_pos=c(480,25))

all_types_rnaseq = all_types
all_types_rnaseq$data = as.matrix(cbind(naive_CD4_PA_genes$RNAseq,
                                        Th1_PA_genes$RNAseq,
                                        Th2_PA_genes$RNAseq,
                                        Th17_PA_genes$RNAseq,
                                        Treg_PA_genes$RNAseq,
                                        Tfollicular_PA_genes$RNAseq,
                                        CD4_CM_PA_genes$RNAseq,
                                        CD4_EM_PA_genes$RNAseq))
all_types_rnaseq$metadata = all_types$metadata[colnames(all_types_rnaseq$data),]
all_types_rnaseq=specmine::convert_to_factor(all_types_rnaseq, 'cell_type')
cluster_analysis(all_types_rnaseq, scale_data=T, color_leafs='cell_type', label_leafs='GSE', leg_pos=c(80,25))
anova_all_types_rnaseq = specmine::aov_all_vars(all_types_rnaseq, 'cell_type')
sig_genes_rnaseq = rownames(anova_all_types_rnaseq)[anova_all_types_rnaseq$fdr<0.05]
heatmap(all_types_rnaseq$data[sig_genes_rnaseq,], Colv=NA, scale='none', main='RNAseq')
heatmap(all_types_rnaseq$data[sig_genes_microarray,], Colv=NA, scale='none', main='RNAseq with Microarray sig genes')

all_types_microarray = all_types
all_types_microarray$data = as.matrix(cbind(naive_CD4_PA_genes$Microarray, Th1_PA_genes$Microarray,Th2_PA_genes$Microarray))
all_types_microarray$metadata = all_types$metadata[colnames(all_types_microarray$data),]
all_types_microarray=specmine::remove_variables_by_nas(all_types_microarray)
all_types_microarray=specmine::convert_to_factor(all_types_microarray, 'cell_type')
anova_all_types_microarray = specmine::aov_all_vars(all_types_microarray, 'cell_type')
sig_genes_microarray = na.omit(rownames(anova_all_types_microarray)[anova_all_types_microarray$fdr<0.05])
heatmap(all_types_microarray$data[sig_genes_microarray,], Colv=NA, scale='none', main='Microarray')
heatmap(all_types_microarray$data[sig_genes_rnaseq,], Colv=NA, scale='none', main='Microarray with RNAseq sig genes')

all_types = specmine::remove_variables_by_nas(all_types)
all_types = specmine::convert_to_factor(all_types, 'cell_type')
anova_all_types = specmine::aov_all_vars(all_types, 'cell_type')
sig_genes = na.omit(rownames(anova_all_types)[anova_all_types$fdr<0.05])
heatmap(all_types$data[sig_genes,], Colv=NA, scale='none', main='Both')


all_types_reactions = cbind(naive_CD4_PA_reactions, Th1_PA_reactions, Th2_PA_reactions)
all_types_reactions_filter = c()
for(reaction in 1:dim(all_types_reactions)[1]){
  if(sum(duplicated(all_types_reactions[reaction,]))==dim(all_types_reactions)[2]-1) all_types_reactions_filter = c(all_types_reactions_filter, F)
  else all_types_reactions_filter = c(all_types_reactions_filter, T)
}
all_types_reactions_filtered = all_types_reactions[all_types_reactions_filter,]
heatmap(na.omit(all_types_reactions_filtered), Colv=NA, scale='none', main='Reactions')

