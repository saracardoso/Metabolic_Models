#Directories:
base_dir = '~/Documents/PhD'

Metabolic_Models_Repo = paste(base_dir, 'Metabolic_Models', sep='/')
Metabolic_Models_Repo_general = paste(Metabolic_Models_Repo, 'general', sep='/')
Metabolic_Models_Repo_general_utils = paste(Metabolic_Models_Repo_general, 'utils', sep='/')
Metabolic_Models_Repo_general_R = paste(Metabolic_Models_Repo_general, 'R', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS = paste(Metabolic_Models_Repo, 'reconstruction_process_T_CELLS', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD8 = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS, 'naive_CD8', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD8_PA_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD8, 'PA_calls', sep='/')

T_Cell_Data_Repo = paste(base_dir, 'T_Cell_Data_Repo', sep='/')
T_Cell_Data_Repo_Data_Files = paste(T_Cell_Data_Repo, 'Data_Files', sep='/')
T_Cell_Data_Repo_Data_Files_naive_CD8 = paste(T_Cell_Data_Repo_Data_Files, 'naive_CD8', sep='/')
T_Cell_Data_Repo_Data_Files_naive_CD8_RNA_seq = paste(T_Cell_Data_Repo_Data_Files_naive_CD8, 'RNA_seq', sep='/')
T_Cell_Data_Repo_Data_Files_naive_CD8_Microarray = paste(T_Cell_Data_Repo_Data_Files_naive_CD8, 'Microarray', sep='/')
T_Cell_Data_Repo_Data_Files_naive_CD8_Proteomics = paste(T_Cell_Data_Repo_Data_Files_naive_CD8, 'Proteomics', sep='/')

#Files:
PA_calls_file = paste(Metabolic_Models_Repo_general_R, 'PA_calls.R', sep='/')
recon3DModel_genes_file = paste(Metabolic_Models_Repo_general_utils, 'recon3DModel_genes.csv', sep='/')
recon3DModel_GPR_file = paste(Metabolic_Models_Repo_general_utils, 'recon3DModel_GPR.txt', sep='/')

naive_CD8_microarray_data_file = paste(T_Cell_Data_Repo_Data_Files_naive_CD8_Microarray, 'naive_CD8_microarray_data.csv', sep='/')
naive_CD8_rnaseq_data_file = paste(T_Cell_Data_Repo_Data_Files_naive_CD8_RNA_seq, 'naive_CD8_rnaseq_data.csv', sep='/')
naive_CD8_proteomics_data_file = paste(T_Cell_Data_Repo_Data_Files_naive_CD8_Proteomics, 'naive_CD8_proteomics_data.csv', sep='/')
naive_CD8_microarray_metadata_file = paste(T_Cell_Data_Repo_Data_Files_naive_CD8_Microarray, 'naive_CD8_microarray_metadata.csv', sep='/')
naive_CD8_rnaseq_metadata_file = paste(T_Cell_Data_Repo_Data_Files_naive_CD8_RNA_seq, 'naive_CD8_rnaseq_metadata.csv', sep='/')

naive_CD8_gene_calls_RNAseq_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD8_PA_calls, 'naive_CD8_gene_calls_RNAseq.csv', sep='/')
naive_CD8_gene_calls_Microarray_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD8_PA_calls, 'naive_CD8_gene_calls_Microarray.csv', sep='/')
naive_CD8_gene_calls_proteomics_ImmProt_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD8_PA_calls, 'naive_CD8_gene_calls_Proteomics_ImmProt.csv', sep='/')
naive_CD8_final_gene_calls_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD8_PA_calls, 'naive_CD8_final_gene_calls.csv', sep='/')
naive_CD8_reaction_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD8_PA_calls, 'naive_CD8_reaction_calls.csv', sep='/')


#Source:
source(PA_calls_file)
recon3DModel_genes = read.csv(recon3DModel_genes_file, stringsAsFactors = F)
recon3D_gpr_rules = read.table(recon3DModel_GPR_file, sep='\t', stringsAsFactors = F)
recon3D_gpr_rules_filtered = recon3D_gpr_rules[!recon3D_gpr_rules[,2]=='',] #Remove reactions without GPRs
recon3D_gpr_rules_filtered[,2] = gsub('[.][1-9]+', '',recon3D_gpr_rules_filtered[,2]) #Remove 'version' .1, .2, ...





#--- Read expression data for PA calls ---

naive_CD8_expression_data = prepare_data_pa_calls(transcriptomics_data_files = list(RNAseq=naive_CD8_rnaseq_data_file, Microarray=naive_CD8_microarray_data_file),
                                               transcriptomics_metadata_files = list(RNAseq=naive_CD8_rnaseq_metadata_file, Microarray=naive_CD8_microarray_metadata_file),
                                               proteomics_data_files = list(ImmProt=naive_CD8_proteomics_data_file),
                                               proteomics_metadata_files = NULL)


#--- Gene Calls ---

naive_CD8_transcriptomics_PA_calls = transcriptomics_PA_calls(naive_CD8_expression_data$Transcriptomics)
View(naive_CD8_transcriptomics_PA_calls$calls_per_omics)
View(naive_CD8_transcriptomics_PA_calls$calls_per_sample$RNAseq)
View(naive_CD8_transcriptomics_PA_calls$calls_per_sample$Microarray)

naive_CD8_proteomics_PA_calls = proteomics_PA_calls(naive_CD8_expression_data$Proteomics)
View(naive_CD8_proteomics_PA_calls$calls_per_omics)
View(naive_CD8_proteomics_PA_calls$calls_per_sample$ImmProt)

naive_CD8_gene_PA_calls = final_gene_PA_calls(naive_CD8_transcriptomics_PA_calls, naive_CD8_proteomics_PA_calls)
View(naive_CD8_gene_PA_calls$final_gene_calls)

write.csv(naive_CD8_transcriptomics_PA_calls$calls_per_sample$RNAseq, naive_CD8_gene_calls_RNAseq_new_file)
write.csv(naive_CD8_transcriptomics_PA_calls$calls_per_sample$Microarray, naive_CD8_gene_calls_Microarray_new_file)
write.csv(naive_CD8_proteomics_PA_calls$calls_per_sample$ImmProt, naive_CD8_gene_calls_proteomics_ImmProt_new_file)
write.csv(naive_CD8_gene_PA_calls$final_gene_calls, naive_CD8_final_gene_calls_new_file)


#-- Reaction Calls --

naive_CD8_PA_reactions = PA_reaction_calls(naive_CD8_gene_PA_calls$final_gene_calls)
View(naive_CD8_PA_reactions)

write.table(naive_CD8_PA_reactions, naive_CD8_reaction_calls, sep=',', col.names=F)


#-- Analysis of PA Calls --

#Heatmap with rnaseq and microarray samples and the final gene calls:
heatmap_gene_calls(naive_CD8_gene_PA_calls, plot_tilte='naive_CD8 | PA Gene Calls',
                   sample_colors_legend_location=c(65,0), cell_colours_legend_location=c(65,10))

naive_CD8_data_to_plot = list()
naive_CD8_data_to_plot$data = cbind(naive_CD8_gene_PA_calls$calls_per_sample$Transcriptomics$RNAseq,
                                 naive_CD8_gene_PA_calls$calls_per_sample$Transcriptomics$Microarray,
                                 naive_CD8_gene_PA_calls$calls_per_sample$Proteomics$ImmProt)
naive_CD8_data_to_plot$metadata = data.frame(GSE=c(naive_CD8_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                                naive_CD8_expression_data$Transcriptomics$Microarray$metadata$GSE,
                                                c(rep('Steady_State', 4),rep('Activated', 3))),
                                          data_type=c(rep('RNAseq', dim(naive_CD8_expression_data$Transcriptomics$RNAseq$data)[2]),
                                                      rep('Microarray', dim(naive_CD8_expression_data$Transcriptomics$Microarray$data)[2]),
                                                      rep('ImmProt', dim(naive_CD8_expression_data$Proteomics$ImmProt$data)[2])),
                                          row.names = c(rownames(naive_CD8_expression_data$Transcriptomics$RNAseq$metadata),
                                                        rownames(naive_CD8_expression_data$Transcriptomics$Microarray$metadata),
                                                        colnames(naive_CD8_gene_PA_calls$calls_per_sample$Proteomics$ImmProt)),
                                          stringsAsFactors = T)
cluster_analysis(naive_CD8_data_to_plot, color_leafs='GSE', label_leafs='data_type', leg_pos=c(25,25)) #Color by GSE dataset and sample names by data type

#How many genes with score 0, 1 and 2:
sum(na.omit(naive_CD8_gene_PA_calls$final_gene_calls$Final==0)) # 1282
sum(na.omit(naive_CD8_gene_PA_calls$final_gene_calls$Final==1)) # 254
sum(na.omit(naive_CD8_gene_PA_calls$final_gene_calls$Final==2)) # 346

#How many reactions with score 0, 1 and 2:
sum(na.omit(naive_CD8_PA_reactions==0)) # 3444
sum(na.omit(naive_CD8_PA_reactions==1)) # 1226
sum(na.omit(naive_CD8_PA_reactions==2)) # 1268

#GO analysis:
# Get the genes from all the reactions present:
naive_CD8_reactions = names(na.omit(naive_CD8_PA_reactions))[na.omit(naive_CD8_PA_reactions)!=0]
naive_CD8_genes_from_reactions=recon3D_gpr_rules_filtered[match(naive_CD8_reactions, recon3D_gpr_rules_filtered[,1]),2]
naive_CD8_genes_from_reactions = unique(unlist(strsplit(unlist(strsplit(gsub('[)]','', gsub('[(]','', na.omit(naive_CD8_genes_from_reactions))), ' or ')), ' and ')))
naive_CD8_genes_from_reactions = naive_CD8_gene_PA_calls$final_gene_calls[match(naive_CD8_genes_from_reactions, rownames(naive_CD8_gene_PA_calls$final_gene_calls)),]
naive_CD8_genes_from_reactions = rownames(naive_CD8_genes_from_reactions)[naive_CD8_genes_from_reactions[,'Final']!=0]
naive_CD8_genes_from_reactions = recon3DModel_genes$ensembl[match(naive_CD8_genes_from_reactions, recon3DModel_genes$entrez)]
# Universe of genes: all genes in recon3D model:
universe_genes = unique(recon3DModel_genes$ensembl)
# GO analysis:
naive_CD8_GO_analysis_reactions = GO_analysis(naive_CD8_genes_from_reactions, universe_genes)
naive_CD8_GO_tables_reactions_sig = significant_GOs(naive_CD8_GO_analysis_reactions)
View(naive_CD8_GO_tables_reactions_sig$BP)
View(naive_CD8_GO_tables_reactions_sig$CC)
View(naive_CD8_GO_tables_reactions_sig$MF)
