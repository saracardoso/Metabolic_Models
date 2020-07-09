#Directories:
base_dir = '~/Documents/PhD'

Metabolic_Models_Repo = paste(base_dir, 'Metabolic_Models', sep='/')
Metabolic_Models_Repo_general = paste(Metabolic_Models_Repo, 'general', sep='/')
Metabolic_Models_Repo_general_utils = paste(Metabolic_Models_Repo_general, 'utils', sep='/')
Metabolic_Models_Repo_general_R = paste(Metabolic_Models_Repo_general, 'R', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS = paste(Metabolic_Models_Repo, 'reconstruction_process_T_CELLS', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_Th1 = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS, 'Th1', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_Th1_PA_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Th1, 'PA_calls', sep='/')

T_Cell_Data_Repo = paste(base_dir, 'T_Cell_Data_Repo', sep='/')
T_Cell_Data_Repo_Data_Files = paste(T_Cell_Data_Repo, 'Data_Files', sep='/')
T_Cell_Data_Repo_Data_Files_Th1 = paste(T_Cell_Data_Repo_Data_Files, 'Th1', sep='/')
T_Cell_Data_Repo_Data_Files_Th1_RNA_seq = paste(T_Cell_Data_Repo_Data_Files_Th1, 'RNA_seq', sep='/')
T_Cell_Data_Repo_Data_Files_Th1_Microarray = paste(T_Cell_Data_Repo_Data_Files_Th1, 'Microarray', sep='/')
T_Cell_Data_Repo_Data_Files_Th1_Proteomics = paste(T_Cell_Data_Repo_Data_Files_Th1, 'Proteomics', sep='/')

#Files:
PA_calls_file = paste(Metabolic_Models_Repo_general_R, 'PA_calls.R', sep='/')
recon3DModel_genes_file = paste(Metabolic_Models_Repo_general_utils, 'recon3DModel_genes.csv', sep='/')
recon3DModel_GPR_file = paste(Metabolic_Models_Repo_general_utils, 'recon3DModel_GPR.txt', sep='/')

Th1_microarray_data_file = paste(T_Cell_Data_Repo_Data_Files_Th1_Microarray, 'Th1_microarray_data.csv', sep='/')
Th1_rnaseq_data_file = paste(T_Cell_Data_Repo_Data_Files_Th1_RNA_seq, 'Th1_rnaseq_data.csv', sep='/')
Th1_proteomics_data_file = paste(T_Cell_Data_Repo_Data_Files_Th1_Proteomics, 'Th1_proteomics_data.csv', sep='/')
Th1_microarray_metadata_file = paste(T_Cell_Data_Repo_Data_Files_Th1_Microarray, 'Th1_microarray_metadata.csv', sep='/')
Th1_rnaseq_metadata_file = paste(T_Cell_Data_Repo_Data_Files_Th1_RNA_seq, 'Th1_rnaseq_metadata.csv', sep='/')

Th1_gene_calls_RNAseq_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Th1_PA_calls, 'Th1_gene_calls_RNAseq.csv', sep='/')
Th1_gene_calls_Microarray_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Th1_PA_calls, 'Th1_gene_calls_Microarray.csv', sep='/')
Th1_gene_calls_proteomics_ImmProt_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Th1_PA_calls, 'Th1_gene_calls_Proteomics_ImmProt.csv', sep='/')
Th1_final_gene_calls_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Th1_PA_calls, 'Th1_final_gene_calls.csv', sep='/')
Th1_reaction_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Th1_PA_calls, 'Th1_reaction_calls.csv', sep='/')


#Source:
source(PA_calls_file)
recon3DModel_genes = read.csv(recon3DModel_genes_file, stringsAsFactors = F)
recon3D_gpr_rules = read.table(recon3DModel_GPR_file, sep='\t', stringsAsFactors = F)
recon3D_gpr_rules_filtered = recon3D_gpr_rules[!recon3D_gpr_rules[,2]=='',] #Remove reactions without GPRs
recon3D_gpr_rules_filtered[,2] = gsub('[.][1-9]+', '',recon3D_gpr_rules_filtered[,2]) #Remove 'version' .1, .2, ...





#--- Read expression data for PA calls ---

Th1_expression_data = prepare_data_pa_calls(transcriptomics_data_files = list(RNAseq=Th1_rnaseq_data_file, Microarray=Th1_microarray_data_file),
                                            transcriptomics_metadata_files = list(RNAseq=Th1_rnaseq_metadata_file, Microarray=Th1_microarray_metadata_file),
                                            proteomics_data_files = list(ImmProt=Th1_proteomics_data_file),
                                            proteomics_metadata_files = NULL)


#--- Gene Calls ---

Th1_transcriptomics_PA_calls = transcriptomics_PA_calls(Th1_expression_data$Transcriptomics)
View(Th1_transcriptomics_PA_calls$calls_per_omics)
View(Th1_transcriptomics_PA_calls$calls_per_sample$RNAseq)
View(Th1_transcriptomics_PA_calls$calls_per_sample$Microarray)

Th1_proteomics_PA_calls = proteomics_PA_calls(Th1_expression_data$Proteomics)
View(Th1_proteomics_PA_calls$calls_per_omics)
View(Th1_proteomics_PA_calls$calls_per_sample$ImmProt)

Th1_gene_PA_calls = final_gene_PA_calls(Th1_transcriptomics_PA_calls, Th1_proteomics_PA_calls)
View(Th1_gene_PA_calls$final_gene_calls)

write.csv(Th1_transcriptomics_PA_calls$calls_per_sample$RNAseq, Th1_gene_calls_RNAseq_new_file)
write.csv(Th1_transcriptomics_PA_calls$calls_per_sample$Microarray, Th1_gene_calls_Microarray_new_file)
write.csv(Th1_proteomics_PA_calls$calls_per_sample$ImmProt, Th1_gene_calls_proteomics_ImmProt_new_file)
write.csv(Th1_gene_PA_calls$final_gene_calls, Th1_final_gene_calls_new_file)


#-- Reaction Calls --

Th1_PA_reactions = PA_reaction_calls(Th1_gene_PA_calls$final_gene_calls)
View(Th1_PA_reactions)

write.table(Th1_PA_reactions, Th1_reaction_calls, sep=',', col.names=F)


#-- Analysis of PA Calls --

#Heatmap with rnaseq and microarray samples and the final gene calls:
heatmap_gene_calls(Th1_gene_PA_calls, plot_tilte='Th1 | PA Gene Calls',
                   sample_colors_legend_location=c(65,0), cell_colours_legend_location=c(65,10))

Th1_data_to_plot = list()
Th1_data_to_plot$data = cbind(Th1_gene_PA_calls$calls_per_sample$Transcriptomics$RNAseq,
                                    Th1_gene_PA_calls$calls_per_sample$Transcriptomics$Microarray,
                                    Th1_gene_PA_calls$calls_per_sample$Proteomics$ImmProt)
Th1_data_to_plot$metadata = data.frame(GSE=c(Th1_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                                   Th1_expression_data$Transcriptomics$Microarray$metadata$GSE,
                                                   c(rep('Steady_State', 4))),
                                             data_type=c(rep('RNAseq', dim(Th1_expression_data$Transcriptomics$RNAseq$data)[2]),
                                                         rep('Microarray', dim(Th1_expression_data$Transcriptomics$Microarray$data)[2]),
                                                         rep('ImmProt', dim(Th1_expression_data$Proteomics$ImmProt$data)[2])),
                                             row.names = c(rownames(Th1_expression_data$Transcriptomics$RNAseq$metadata),
                                                           rownames(Th1_expression_data$Transcriptomics$Microarray$metadata),
                                                           colnames(Th1_gene_PA_calls$calls_per_sample$Proteomics$ImmProt)),
                                             stringsAsFactors = T)
cluster_analysis(Th1_data_to_plot, color_leafs='GSE', label_leafs='data_type', leg_pos=c(25,25)) #Color by GSE dataset and sample names by data type

#How many genes with score 0, 1 and 2:
sum(na.omit(Th1_gene_PA_calls$final_gene_calls$Final==0)) # 986
sum(na.omit(Th1_gene_PA_calls$final_gene_calls$Final==1)) # 428
sum(na.omit(Th1_gene_PA_calls$final_gene_calls$Final==2)) # 468

#How many reactions with score 0, 1 and 2:
sum(na.omit(Th1_PA_reactions==0)) # 2435
sum(na.omit(Th1_PA_reactions==1)) # 1346
sum(na.omit(Th1_PA_reactions==2)) # 2157

#GO analysis:
# Get the genes from all the reactions present:
Th1_reactions = names(na.omit(Th1_PA_reactions))[na.omit(Th1_PA_reactions)!=0]
Th1_genes_from_reactions=recon3D_gpr_rules_filtered[match(Th1_reactions, recon3D_gpr_rules_filtered[,1]),2]
Th1_genes_from_reactions = unique(unlist(strsplit(unlist(strsplit(gsub('[)]','', gsub('[(]','', na.omit(Th1_genes_from_reactions))), ' or ')), ' and ')))
Th1_genes_from_reactions = Th1_gene_PA_calls$final_gene_calls[match(Th1_genes_from_reactions, rownames(Th1_gene_PA_calls$final_gene_calls)),]
Th1_genes_from_reactions = rownames(Th1_genes_from_reactions)[Th1_genes_from_reactions[,'Final']!=0]
Th1_genes_from_reactions = recon3DModel_genes$ensembl[match(Th1_genes_from_reactions, recon3DModel_genes$entrez)]
# Universe of genes: all genes in recon3D model:
universe_genes = unique(recon3DModel_genes$ensembl)
# GO analysis:
Th1_GO_analysis_reactions = GO_analysis(Th1_genes_from_reactions, universe_genes)
Th1_GO_tables_reactions_sig = significant_GOs(Th1_GO_analysis_reactions)
View(Th1_GO_tables_reactions_sig$BP)
View(Th1_GO_tables_reactions_sig$CC)
View(Th1_GO_tables_reactions_sig$MF)
