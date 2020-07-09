#Directories:
base_dir = '~/Documents/PhD'

Metabolic_Models_Repo = paste(base_dir, 'Metabolic_Models', sep='/')
Metabolic_Models_Repo_general = paste(Metabolic_Models_Repo, 'general', sep='/')
Metabolic_Models_Repo_general_utils = paste(Metabolic_Models_Repo_general, 'utils', sep='/')
Metabolic_Models_Repo_general_R = paste(Metabolic_Models_Repo_general, 'R', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS = paste(Metabolic_Models_Repo, 'reconstruction_process_T_CELLS', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_CD4_memory = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS, 'CD4_memory', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_CD4_memory_PA_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_CD4_memory, 'PA_calls', sep='/')

T_Cell_Data_Repo = paste(base_dir, 'T_Cell_Data_Repo', sep='/')
T_Cell_Data_Repo_Data_Files = paste(T_Cell_Data_Repo, 'Data_Files', sep='/')
T_Cell_Data_Repo_Data_Files_CD4_memory = paste(T_Cell_Data_Repo_Data_Files, 'CD4_memory', sep='/')
T_Cell_Data_Repo_Data_Files_CD4_memory_Microarray = paste(T_Cell_Data_Repo_Data_Files_CD4_memory, 'Microarray', sep='/')
T_Cell_Data_Repo_Data_Files_CD4_memory_Proteomics = paste(T_Cell_Data_Repo_Data_Files_CD4_memory, 'Proteomics', sep='/')

#Files:
PA_calls_file = paste(Metabolic_Models_Repo_general_R, 'PA_calls.R', sep='/')
recon3DModel_genes_file = paste(Metabolic_Models_Repo_general_utils, 'recon3DModel_genes.csv', sep='/')
recon3DModel_GPR_file = paste(Metabolic_Models_Repo_general_utils, 'recon3DModel_GPR.txt', sep='/')

CD4_memory_microarray_data_file = paste(T_Cell_Data_Repo_Data_Files_CD4_memory_Microarray, 'CD4_memory_microarray_data.csv', sep='/')
CD4_memory_microarray_metadata_file = paste(T_Cell_Data_Repo_Data_Files_CD4_memory_Microarray, 'CD4_memory_microarray_metadata.csv', sep='/')

CD4_memory_gene_calls_Microarray_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_CD4_memory_PA_calls, 'CD4_memory_gene_calls_Microarray.csv', sep='/')
CD4_memory_final_gene_calls_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_CD4_memory_PA_calls, 'CD4_memory_final_gene_calls.csv', sep='/')
CD4_memory_reaction_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_CD4_memory_PA_calls, 'CD4_memory_reaction_calls.csv', sep='/')


#Source:
source(PA_calls_file)
recon3DModel_genes = read.csv(recon3DModel_genes_file, stringsAsFactors = F)
recon3D_gpr_rules = read.table(recon3DModel_GPR_file, sep='\t', stringsAsFactors = F)
recon3D_gpr_rules_filtered = recon3D_gpr_rules[!recon3D_gpr_rules[,2]=='',] #Remove reactions without GPRs
recon3D_gpr_rules_filtered[,2] = gsub('[.][1-9]+', '',recon3D_gpr_rules_filtered[,2]) #Remove 'version' .1, .2, ...





#--- Read expression data for PA calls ---

CD4_memory_expression_data = prepare_data_pa_calls(transcriptomics_data_files = list(Microarray=CD4_memory_microarray_data_file),
                                                    transcriptomics_metadata_files = list(Microarray=CD4_memory_microarray_metadata_file),
                                                    proteomics_data_files = NULL,
                                                    proteomics_metadata_files = NULL)


#--- Gene Calls ---

CD4_memory_transcriptomics_PA_calls = transcriptomics_PA_calls(CD4_memory_expression_data$Transcriptomics)
View(CD4_memory_transcriptomics_PA_calls$calls_per_omics)
View(CD4_memory_transcriptomics_PA_calls$calls_per_sample$Microarray)

CD4_memory_gene_PA_calls = final_gene_PA_calls(CD4_memory_transcriptomics_PA_calls, NULL)
View(CD4_memory_gene_PA_calls$final_gene_calls)

write.csv(CD4_memory_transcriptomics_PA_calls$calls_per_sample$Microarray, CD4_memory_gene_calls_Microarray_new_file)
write.csv(CD4_memory_gene_PA_calls$final_gene_calls, CD4_memory_final_gene_calls_new_file)


#-- Reaction Calls --

CD4_memory_PA_reactions = PA_reaction_calls(CD4_memory_gene_PA_calls$final_gene_calls)
View(CD4_memory_PA_reactions)

write.table(CD4_memory_PA_reactions, CD4_memory_reaction_calls, sep=',', col.names=F)


#-- Analysis of PA Calls --

#Heatmap with rnaseq and microarray samples and the final gene calls:
heatmap_gene_calls(CD4_memory_gene_PA_calls, plot_tilte='CD4_memory | PA Gene Calls',
                   sample_colors_legend_location=c(65,0), cell_colours_legend_location=c(65,10))

CD4_memory_data_to_plot = list()
CD4_memory_data_to_plot$data = cbind(CD4_memory_gene_PA_calls$calls_per_sample$Transcriptomics$Microarray)
CD4_memory_data_to_plot$metadata = data.frame(GSE=c(CD4_memory_expression_data$Transcriptomics$Microarray$metadata$GSE),
                                               data_type=c(rep('Microarray', dim(CD4_memory_expression_data$Transcriptomics$Microarray$data)[2])),
                                               row.names = c(rownames(CD4_memory_expression_data$Transcriptomics$Microarray$metadata)),
                                               stringsAsFactors = T)
cluster_analysis(CD4_memory_data_to_plot, color_leafs='GSE', label_leafs='data_type', leg_pos=c(25,25)) #Color by GSE dataset and sample names by data type

#How many genes with score 0, 1 and 2:
sum(na.omit(CD4_memory_gene_PA_calls$final_gene_calls$Final==0)) # 879
sum(na.omit(CD4_memory_gene_PA_calls$final_gene_calls$Final==1)) # 1003
sum(na.omit(CD4_memory_gene_PA_calls$final_gene_calls$Final==2)) # 0

#How many reactions with score 0, 1 and 2:
sum(na.omit(CD4_memory_PA_reactions==0)) # 2122
sum(na.omit(CD4_memory_PA_reactions==1)) # 3816
sum(na.omit(CD4_memory_PA_reactions==2)) # 0

#GO analysis:
# Get the genes from all the reactions present:
CD4_memory_reactions = names(na.omit(CD4_memory_PA_reactions))[na.omit(CD4_memory_PA_reactions)!=0]
CD4_memory_genes_from_reactions=recon3D_gpr_rules_filtered[match(CD4_memory_reactions, recon3D_gpr_rules_filtered[,1]),2]
CD4_memory_genes_from_reactions = unique(unlist(strsplit(unlist(strsplit(gsub('[)]','', gsub('[(]','', na.omit(CD4_memory_genes_from_reactions))), ' or ')), ' and ')))
CD4_memory_genes_from_reactions = CD4_memory_gene_PA_calls$final_gene_calls[match(CD4_memory_genes_from_reactions, rownames(CD4_memory_gene_PA_calls$final_gene_calls)),]
CD4_memory_genes_from_reactions = rownames(CD4_memory_genes_from_reactions)[CD4_memory_genes_from_reactions[,'Final']!=0]
CD4_memory_genes_from_reactions = recon3DModel_genes$ensembl[match(CD4_memory_genes_from_reactions, recon3DModel_genes$entrez)]
# Universe of genes: all genes in recon3D model:
universe_genes = unique(recon3DModel_genes$ensembl)
# GO analysis:
CD4_memory_GO_analysis_reactions = GO_analysis(CD4_memory_genes_from_reactions, universe_genes)
CD4_memory_GO_tables_reactions_sig = significant_GOs(CD4_memory_GO_analysis_reactions)
View(CD4_memory_GO_tables_reactions_sig$BP)
View(CD4_memory_GO_tables_reactions_sig$CC)
View(CD4_memory_GO_tables_reactions_sig$MF)
