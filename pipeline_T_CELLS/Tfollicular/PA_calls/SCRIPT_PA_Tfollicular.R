#Directories:
base_dir = '~/Documents/PhD'

Metabolic_Models_Repo = paste(base_dir, 'Metabolic_Models', sep='/')
Metabolic_Models_Repo_general = paste(Metabolic_Models_Repo, 'general', sep='/')
Metabolic_Models_Repo_general_utils = paste(Metabolic_Models_Repo_general, 'utils', sep='/')
Metabolic_Models_Repo_general_R = paste(Metabolic_Models_Repo_general, 'R', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS = paste(Metabolic_Models_Repo, 'reconstruction_process_T_CELLS', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_Tfollicular = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS, 'Tfollicular', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_Tfollicular_PA_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Tfollicular, 'PA_calls', sep='/')

T_Cell_Data_Repo = paste(base_dir, 'T_Cell_Data_Repo', sep='/')
T_Cell_Data_Repo_Data_Files = paste(T_Cell_Data_Repo, 'Data_Files', sep='/')
T_Cell_Data_Repo_Data_Files_Tfollicular = paste(T_Cell_Data_Repo_Data_Files, 'Tfollicular', sep='/')
T_Cell_Data_Repo_Data_Files_Tfollicular_RNA_seq = paste(T_Cell_Data_Repo_Data_Files_Tfollicular, 'RNA_seq', sep='/')
T_Cell_Data_Repo_Data_Files_Tfollicular_Microarray = paste(T_Cell_Data_Repo_Data_Files_Tfollicular, 'Microarray', sep='/')
T_Cell_Data_Repo_Data_Files_Tfollicular_Proteomics = paste(T_Cell_Data_Repo_Data_Files_Tfollicular, 'Proteomics', sep='/')

#Files:
PA_calls_file = paste(Metabolic_Models_Repo_general_R, 'PA_calls.R', sep='/')
recon3DModel_genes_file = paste(Metabolic_Models_Repo_general_utils, 'recon3DModel_genes.csv', sep='/')
recon3DModel_GPR_file = paste(Metabolic_Models_Repo_general_utils, 'recon3DModel_GPR.txt', sep='/')

Tfollicular_microarray_data_file = paste(T_Cell_Data_Repo_Data_Files_Tfollicular_Microarray, 'Tfollicular_microarray_data.csv', sep='/')
Tfollicular_rnaseq_data_file = paste(T_Cell_Data_Repo_Data_Files_Tfollicular_RNA_seq, 'Tfollicular_rnaseq_data.csv', sep='/')
Tfollicular_microarray_metadata_file = paste(T_Cell_Data_Repo_Data_Files_Tfollicular_Microarray, 'Tfollicular_microarray_metadata.csv', sep='/')
Tfollicular_rnaseq_metadata_file = paste(T_Cell_Data_Repo_Data_Files_Tfollicular_RNA_seq, 'Tfollicular_rnaseq_metadata.csv', sep='/')

Tfollicular_gene_calls_RNAseq_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Tfollicular_PA_calls, 'Tfollicular_gene_calls_RNAseq.csv', sep='/')
Tfollicular_gene_calls_Microarray_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Tfollicular_PA_calls, 'Tfollicular_gene_calls_Microarray.csv', sep='/')
Tfollicular_final_gene_calls_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Tfollicular_PA_calls, 'Tfollicular_final_gene_calls.csv', sep='/')
Tfollicular_reaction_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Tfollicular_PA_calls, 'Tfollicular_reaction_calls.csv', sep='/')


#Source:
source(PA_calls_file)
recon3DModel_genes = read.csv(recon3DModel_genes_file, stringsAsFactors = F)
recon3D_gpr_rules = read.table(recon3DModel_GPR_file, sep='\t', stringsAsFactors = F)
recon3D_gpr_rules_filtered = recon3D_gpr_rules[!recon3D_gpr_rules[,2]=='',] #Remove reactions without GPRs
recon3D_gpr_rules_filtered[,2] = gsub('[.][1-9]+', '',recon3D_gpr_rules_filtered[,2]) #Remove 'version' .1, .2, ...





#--- Read expression data for PA calls ---

Tfollicular_expression_data = prepare_data_pa_calls(transcriptomics_data_files = list(RNAseq=Tfollicular_rnaseq_data_file, Microarray=Tfollicular_microarray_data_file),
                                             transcriptomics_metadata_files = list(RNAseq=Tfollicular_rnaseq_metadata_file, Microarray=Tfollicular_microarray_metadata_file),
                                             proteomics_data_files = NULL,
                                             proteomics_metadata_files = NULL)


#--- Gene Calls ---

Tfollicular_transcriptomics_PA_calls = transcriptomics_PA_calls(Tfollicular_expression_data$Transcriptomics)
View(Tfollicular_transcriptomics_PA_calls$calls_per_omics)
View(Tfollicular_transcriptomics_PA_calls$calls_per_sample$RNAseq)
View(Tfollicular_transcriptomics_PA_calls$calls_per_sample$Microarray)

Tfollicular_gene_PA_calls = final_gene_PA_calls(Tfollicular_transcriptomics_PA_calls, NULL)
View(Tfollicular_gene_PA_calls$final_gene_calls)

write.csv(Tfollicular_transcriptomics_PA_calls$calls_per_sample$RNAseq, Tfollicular_gene_calls_RNAseq_new_file)
write.csv(Tfollicular_transcriptomics_PA_calls$calls_per_sample$Microarray, Tfollicular_gene_calls_Microarray_new_file)
write.csv(Tfollicular_gene_PA_calls$final_gene_calls, Tfollicular_final_gene_calls_new_file)


#-- Reaction Calls --

Tfollicular_PA_reactions = PA_reaction_calls(Tfollicular_gene_PA_calls$final_gene_calls)
View(Tfollicular_PA_reactions)

write.table(Tfollicular_PA_reactions, Tfollicular_reaction_calls, sep=',', col.names=F)


#-- Analysis of PA Calls --

#Heatmap with rnaseq and microarray samples and the final gene calls:
heatmap_gene_calls(Tfollicular_gene_PA_calls, plot_tilte='Tfollicular | PA Gene Calls',
                   sample_colors_legend_location=c(65,0), cell_colours_legend_location=c(65,10))

Tfollicular_data_to_plot = list()
Tfollicular_data_to_plot$data = cbind(Tfollicular_gene_PA_calls$calls_per_sample$Transcriptomics$RNAseq,
                               Tfollicular_gene_PA_calls$calls_per_sample$Transcriptomics$Microarray)
Tfollicular_data_to_plot$metadata = data.frame(GSE=c(Tfollicular_expression_data$Transcriptomics$RNAseq$metadata$GSE,
                                              Tfollicular_expression_data$Transcriptomics$Microarray$metadata$GSE),
                                        data_type=c(rep('RNAseq', dim(Tfollicular_expression_data$Transcriptomics$RNAseq$data)[2]),
                                                    rep('Microarray', dim(Tfollicular_expression_data$Transcriptomics$Microarray$data)[2])),
                                        row.names = c(rownames(Tfollicular_expression_data$Transcriptomics$RNAseq$metadata),
                                                      rownames(Tfollicular_expression_data$Transcriptomics$Microarray$metadata),
                                                      colnames(Tfollicular_gene_PA_calls$calls_per_sample$Proteomics$ImmProt)),
                                        stringsAsFactors = T)
cluster_analysis(Tfollicular_data_to_plot, color_leafs='GSE', label_leafs='data_type', leg_pos=c(25,25)) #Color by GSE dataset and sample names by data type

#How many genes with score 0, 1 and 2:
sum(na.omit(Tfollicular_gene_PA_calls$final_gene_calls$Final==0)) # 1050
sum(na.omit(Tfollicular_gene_PA_calls$final_gene_calls$Final==1)) # 540
sum(na.omit(Tfollicular_gene_PA_calls$final_gene_calls$Final==2)) # 292

#How many reactions with score 0, 1 and 2:
sum(na.omit(Tfollicular_PA_reactions==0)) # 2552
sum(na.omit(Tfollicular_PA_reactions==1)) # 2242
sum(na.omit(Tfollicular_PA_reactions==2)) # 1144

#GO analysis:
# Get the genes from all the reactions present:
Tfollicular_reactions = names(na.omit(Tfollicular_PA_reactions))[na.omit(Tfollicular_PA_reactions)!=0]
Tfollicular_genes_from_reactions=recon3D_gpr_rules_filtered[match(Tfollicular_reactions, recon3D_gpr_rules_filtered[,1]),2]
Tfollicular_genes_from_reactions = unique(unlist(strsplit(unlist(strsplit(gsub('[)]','', gsub('[(]','', na.omit(Tfollicular_genes_from_reactions))), ' or ')), ' and ')))
Tfollicular_genes_from_reactions = Tfollicular_gene_PA_calls$final_gene_calls[match(Tfollicular_genes_from_reactions, rownames(Tfollicular_gene_PA_calls$final_gene_calls)),]
Tfollicular_genes_from_reactions = rownames(Tfollicular_genes_from_reactions)[Tfollicular_genes_from_reactions[,'Final']!=0]
Tfollicular_genes_from_reactions = recon3DModel_genes$ensembl[match(Tfollicular_genes_from_reactions, recon3DModel_genes$entrez)]
# Universe of genes: all genes in recon3D model:
universe_genes = unique(recon3DModel_genes$ensembl)
# GO analysis:
Tfollicular_GO_analysis_reactions = GO_analysis(Tfollicular_genes_from_reactions, universe_genes)
Tfollicular_GO_tables_reactions_sig = significant_GOs(Tfollicular_GO_analysis_reactions)
View(Tfollicular_GO_tables_reactions_sig$BP)
View(Tfollicular_GO_tables_reactions_sig$CC)
View(Tfollicular_GO_tables_reactions_sig$MF)
