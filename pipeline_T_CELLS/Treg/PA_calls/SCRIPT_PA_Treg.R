#Directories:
base_dir = '~/Documents/PhD'

Metabolic_Models_Repo = paste(base_dir, 'Metabolic_Models', sep='/')
Metabolic_Models_Repo_general = paste(Metabolic_Models_Repo, 'general', sep='/')
Metabolic_Models_Repo_general_utils = paste(Metabolic_Models_Repo_general, 'utils', sep='/')
Metabolic_Models_Repo_general_R = paste(Metabolic_Models_Repo_general, 'R', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS = paste(Metabolic_Models_Repo, 'reconstruction_process_T_CELLS', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_Treg = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS, 'Treg', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_Treg_PA_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Treg, 'PA_calls', sep='/')

T_Cell_Data_Repo = paste(base_dir, 'T_Cell_Data_Repo', sep='/')
T_Cell_Data_Repo_Data_Files = paste(T_Cell_Data_Repo, 'Data_Files', sep='/')
T_Cell_Data_Repo_Data_Files_Treg = paste(T_Cell_Data_Repo_Data_Files, 'Treg', sep='/')
T_Cell_Data_Repo_Data_Files_Treg_RNA_seq = paste(T_Cell_Data_Repo_Data_Files_Treg, 'RNA_seq', sep='/')
T_Cell_Data_Repo_Data_Files_Treg_Microarray = paste(T_Cell_Data_Repo_Data_Files_Treg, 'Microarray', sep='/')
T_Cell_Data_Repo_Data_Files_Treg_Proteomics = paste(T_Cell_Data_Repo_Data_Files_Treg, 'Proteomics', sep='/')

#Files:
PA_calls_file = paste(Metabolic_Models_Repo_general_R, 'PA_calls.R', sep='/')
recon3DModel_genes_file = paste(Metabolic_Models_Repo_general_utils, 'recon3DModel_genes.csv', sep='/')
recon3DModel_GPR_file = paste(Metabolic_Models_Repo_general_utils, 'recon3DModel_GPR.txt', sep='/')

Treg_microarray_data_file = paste(T_Cell_Data_Repo_Data_Files_Treg_Microarray, 'Treg_microarray_data.csv', sep='/')
Treg_rnaseq_data_file = paste(T_Cell_Data_Repo_Data_Files_Treg_RNA_seq, 'Treg_rnaseq_data.csv', sep='/')
Treg_proteomics_data_file = paste(T_Cell_Data_Repo_Data_Files_Treg_Proteomics, 'Treg_proteomics_data.csv', sep='/')
Treg_microarray_metadata_file = paste(T_Cell_Data_Repo_Data_Files_Treg_Microarray, 'Treg_microarray_metadata.csv', sep='/')
Treg_rnaseq_metadata_file = paste(T_Cell_Data_Repo_Data_Files_Treg_RNA_seq, 'Treg_rnaseq_metadata.csv', sep='/')

Treg_gene_calls_RNAseq_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Treg_PA_calls, 'Treg_gene_calls_RNAseq.csv', sep='/')
Treg_gene_calls_Microarray_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Treg_PA_calls, 'Treg_gene_calls_Microarray.csv', sep='/')
Treg_gene_calls_proteomics_ImmProt_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Treg_PA_calls, 'Treg_gene_calls_Proteomics_ImmProt.csv', sep='/')
Treg_final_gene_calls_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Treg_PA_calls, 'Treg_final_gene_calls.csv', sep='/')
Treg_reaction_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Treg_PA_calls, 'Treg_reaction_calls.csv', sep='/')


#Source:
source(PA_calls_file)
recon3DModel_genes = read.csv(recon3DModel_genes_file, stringsAsFactors = F)
recon3D_gpr_rules = read.table(recon3DModel_GPR_file, sep='\t', stringsAsFactors = F)
recon3D_gpr_rules_filtered = recon3D_gpr_rules[!recon3D_gpr_rules[,2]=='',] #Remove reactions without GPRs
recon3D_gpr_rules_filtered[,2] = gsub('[.][1-9]+', '',recon3D_gpr_rules_filtered[,2]) #Remove 'version' .1, .2, ...





#--- Read expression data for PA calls ---

Treg_expression_data = prepare_data_pa_calls(transcriptomics_data_files = list(RNAseq=Treg_rnaseq_data_file, Microarray=Treg_microarray_data_file),
                                             transcriptomics_metadata_files = list(RNAseq=Treg_rnaseq_metadata_file, Microarray=Treg_microarray_metadata_file),
                                             proteomics_data_files = list(ImmProt=Treg_proteomics_data_file),
                                             proteomics_metadata_files = NULL)


#--- Gene Calls ---

Treg_transcriptomics_PA_calls = transcriptomics_PA_calls(Treg_expression_data$Transcriptomics)
View(Treg_transcriptomics_PA_calls$calls_per_omics)
View(Treg_transcriptomics_PA_calls$calls_per_sample$RNAseq)
View(Treg_transcriptomics_PA_calls$calls_per_sample$Microarray)

Treg_proteomics_PA_calls = proteomics_PA_calls(Treg_expression_data$Proteomics)
View(Treg_proteomics_PA_calls$calls_per_omics)
View(Treg_proteomics_PA_calls$calls_per_sample$ImmProt)

Treg_gene_PA_calls = final_gene_PA_calls(Treg_transcriptomics_PA_calls, Treg_proteomics_PA_calls)
View(Treg_gene_PA_calls$final_gene_calls)

write.csv(Treg_transcriptomics_PA_calls$calls_per_sample$RNAseq, Treg_gene_calls_RNAseq_new_file)
write.csv(Treg_transcriptomics_PA_calls$calls_per_sample$Microarray, Treg_gene_calls_Microarray_new_file)
write.csv(Treg_proteomics_PA_calls$calls_per_sample$ImmProt, Treg_gene_calls_proteomics_ImmProt_new_file)
write.csv(Treg_gene_PA_calls$final_gene_calls, Treg_final_gene_calls_new_file)


#-- Reaction Calls --

Treg_PA_reactions = PA_reaction_calls(Treg_gene_PA_calls$final_gene_calls)
View(Treg_PA_reactions)

write.table(Treg_PA_reactions, Treg_reaction_calls, sep=',', col.names=F)


#-- Analysis of PA Calls --

#Heatmap with rnaseq and microarray samples and the final gene calls:
heatmap_gene_calls(Treg_gene_PA_calls, plot_tilte='Treg | PA Gene Calls',
                   sample_colors_legend_location=c(65,0), cell_colours_legend_location=c(65,10))

Treg_data_to_plot = list()
Treg_data_to_plot$data = cbind(Treg_gene_PA_calls$calls_per_sample$Transcriptomics$RNAseq,
                               Treg_gene_PA_calls$calls_per_sample$Transcriptomics$Microarray,
                               Treg_gene_PA_calls$calls_per_sample$Proteomics$ImmProt)
Treg_data_to_plot$metadata = data.frame(GSE=c(Treg_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                              Treg_expression_data$Transcriptomics$Microarray$metadata$GSE,
                                              c(rep('Steady_State_m', 4),rep('Activated_m', 3),rep('Steady_State_n', 4),rep('Activated_n', 3))),
                                        data_type=c(rep('RNAseq', dim(Treg_expression_data$Transcriptomics$RNAseq$data)[2]),
                                                    rep('Microarray', dim(Treg_expression_data$Transcriptomics$Microarray$data)[2]),
                                                    rep('ImmProt', dim(Treg_expression_data$Proteomics$ImmProt$data)[2])),
                                        row.names = c(rownames(Treg_expression_data$Transcriptomics$RNAseq$metadata),
                                                      rownames(Treg_expression_data$Transcriptomics$Microarray$metadata),
                                                      colnames(Treg_gene_PA_calls$calls_per_sample$Proteomics$ImmProt)),
                                        stringsAsFactors = T)
cluster_analysis(Treg_data_to_plot, color_leafs='GSE', label_leafs='data_type', leg_pos=c(25,25)) #Color by GSE dataset and sample names by data type

#How many genes with score 0, 1 and 2:
sum(na.omit(Treg_gene_PA_calls$final_gene_calls$Final==0)) # 1282
sum(na.omit(Treg_gene_PA_calls$final_gene_calls$Final==1)) # 259
sum(na.omit(Treg_gene_PA_calls$final_gene_calls$Final==2)) # 341

#How many reactions with score 0, 1 and 2:
sum(na.omit(Treg_PA_reactions==0)) # 3575
sum(na.omit(Treg_PA_reactions==1)) # 693
sum(na.omit(Treg_PA_reactions==2)) # 1670

#GO analysis:
# Get the genes from all the reactions present:
Treg_reactions = names(na.omit(Treg_PA_reactions))[na.omit(Treg_PA_reactions)!=0]
Treg_genes_from_reactions=recon3D_gpr_rules_filtered[match(Treg_reactions, recon3D_gpr_rules_filtered[,1]),2]
Treg_genes_from_reactions = unique(unlist(strsplit(unlist(strsplit(gsub('[)]','', gsub('[(]','', na.omit(Treg_genes_from_reactions))), ' or ')), ' and ')))
Treg_genes_from_reactions = Treg_gene_PA_calls$final_gene_calls[match(Treg_genes_from_reactions, rownames(Treg_gene_PA_calls$final_gene_calls)),]
Treg_genes_from_reactions = rownames(Treg_genes_from_reactions)[Treg_genes_from_reactions[,'Final']!=0]
Treg_genes_from_reactions = recon3DModel_genes$ensembl[match(Treg_genes_from_reactions, recon3DModel_genes$entrez)]
# Universe of genes: all genes in recon3D model:
universe_genes = unique(recon3DModel_genes$ensembl)
# GO analysis:
Treg_GO_analysis_reactions = GO_analysis(Treg_genes_from_reactions, universe_genes)
Treg_GO_tables_reactions_sig = significant_GOs(Treg_GO_analysis_reactions)
View(Treg_GO_tables_reactions_sig$BP)
View(Treg_GO_tables_reactions_sig$CC)
View(Treg_GO_tables_reactions_sig$MF)
