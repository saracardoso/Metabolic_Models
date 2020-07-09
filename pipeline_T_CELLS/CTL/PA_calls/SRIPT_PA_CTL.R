#Directories:
base_dir = '~/Documents/PhD'

Metabolic_Models_Repo = paste(base_dir, 'Metabolic_Models', sep='/')
Metabolic_Models_Repo_general = paste(Metabolic_Models_Repo, 'general', sep='/')
Metabolic_Models_Repo_general_utils = paste(Metabolic_Models_Repo_general, 'utils', sep='/')
Metabolic_Models_Repo_general_R = paste(Metabolic_Models_Repo_general, 'R', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS = paste(Metabolic_Models_Repo, 'reconstruction_process_T_CELLS', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_CTL = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS, 'CTL', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_CTL_PA_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_CTL, 'PA_calls', sep='/')

T_Cell_Data_Repo = paste(base_dir, 'T_Cell_Data_Repo', sep='/')
T_Cell_Data_Repo_Data_Files = paste(T_Cell_Data_Repo, 'Data_Files', sep='/')
T_Cell_Data_Repo_Data_Files_CTL = paste(T_Cell_Data_Repo_Data_Files, 'CTL', sep='/')
T_Cell_Data_Repo_Data_Files_CTL_Microarray = paste(T_Cell_Data_Repo_Data_Files_CTL, 'Microarray', sep='/')
T_Cell_Data_Repo_Data_Files_CTL_Proteomics = paste(T_Cell_Data_Repo_Data_Files_CTL, 'Proteomics', sep='/')

#Files:
PA_calls_file = paste(Metabolic_Models_Repo_general_R, 'PA_calls.R', sep='/')
recon3DModel_genes_file = paste(Metabolic_Models_Repo_general_utils, 'recon3DModel_genes.csv', sep='/')
recon3DModel_GPR_file = paste(Metabolic_Models_Repo_general_utils, 'recon3DModel_GPR.txt', sep='/')

CTL_microarray_data_file = paste(T_Cell_Data_Repo_Data_Files_CTL_Microarray, 'CTL_microarray_data.csv', sep='/')
CTL_microarray_metadata_file = paste(T_Cell_Data_Repo_Data_Files_CTL_Microarray, 'CTL_microarray_metadata.csv', sep='/')

CTL_gene_calls_Microarray_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_CTL_PA_calls, 'CTL_gene_calls_Microarray.csv', sep='/')
CTL_final_gene_calls_new_file = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_CTL_PA_calls, 'CTL_final_gene_calls.csv', sep='/')
CTL_reaction_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_CTL_PA_calls, 'CTL_reaction_calls.csv', sep='/')


#Source:
source(PA_calls_file)
recon3DModel_genes = read.csv(recon3DModel_genes_file, stringsAsFactors = F)
recon3D_gpr_rules = read.table(recon3DModel_GPR_file, sep='\t', stringsAsFactors = F)
recon3D_gpr_rules_filtered = recon3D_gpr_rules[!recon3D_gpr_rules[,2]=='',] #Remove reactions without GPRs
recon3D_gpr_rules_filtered[,2] = gsub('[.][1-9]+', '',recon3D_gpr_rules_filtered[,2]) #Remove 'version' .1, .2, ...





#--- Read expression data for PA calls ---

CTL_expression_data = prepare_data_pa_calls(transcriptomics_data_files = list(Microarray=CTL_microarray_data_file),
                                                   transcriptomics_metadata_files = list(Microarray=CTL_microarray_metadata_file),
                                                   proteomics_data_files = NULL,
                                                   proteomics_metadata_files = NULL)


#--- Gene Calls ---

CTL_transcriptomics_PA_calls = transcriptomics_PA_calls(CTL_expression_data$Transcriptomics)
View(CTL_transcriptomics_PA_calls$calls_per_omics)
View(CTL_transcriptomics_PA_calls$calls_per_sample$Microarray)

CTL_gene_PA_calls = final_gene_PA_calls(CTL_transcriptomics_PA_calls, NULL)
View(CTL_gene_PA_calls$final_gene_calls)

write.csv(CTL_transcriptomics_PA_calls$calls_per_sample$Microarray, CTL_gene_calls_Microarray_new_file)
write.csv(CTL_gene_PA_calls$final_gene_calls, CTL_final_gene_calls_new_file)


#-- Reaction Calls --

CTL_PA_reactions = PA_reaction_calls(CTL_gene_PA_calls$final_gene_calls)
View(CTL_PA_reactions)

write.table(CTL_PA_reactions, CTL_reaction_calls, sep=',', col.names=F)


#-- Analysis of PA Calls --

#Heatmap with rnaseq and microarray samples and the final gene calls:
heatmap_gene_calls(CTL_gene_PA_calls, plot_tilte='CTL | PA Gene Calls',
                   sample_colors_legend_location=c(65,0), cell_colours_legend_location=c(65,10))

CTL_data_to_plot = list()
CTL_data_to_plot$data = cbind(CTL_gene_PA_calls$calls_per_sample$Transcriptomics$Microarray)
CTL_data_to_plot$metadata = data.frame(GSE=c(CTL_expression_data$Transcriptomics$Microarray$metadata$GSE),
                                              data_type=c(rep('Microarray', dim(CTL_expression_data$Transcriptomics$Microarray$data)[2])),
                                              row.names = c(rownames(CTL_expression_data$Transcriptomics$Microarray$metadata)),
                                              stringsAsFactors = T)
cluster_analysis(CTL_data_to_plot, color_leafs='GSE', label_leafs='data_type', leg_pos=c(25,25)) #Color by GSE dataset and sample names by data type

#How many genes with score 0, 1 and 2:
sum(na.omit(CTL_gene_PA_calls$final_gene_calls$Final==0)) # 1054
sum(na.omit(CTL_gene_PA_calls$final_gene_calls$Final==1)) # 828
sum(na.omit(CTL_gene_PA_calls$final_gene_calls$Final==2)) # 0

#How many reactions with score 0, 1 and 2:
sum(na.omit(CTL_PA_reactions==0)) # 2749
sum(na.omit(CTL_PA_reactions==1)) # 3189
sum(na.omit(CTL_PA_reactions==2)) # 0

#GO analysis:
# Get the genes from all the reactions present:
CTL_reactions = names(na.omit(CTL_PA_reactions))[na.omit(CTL_PA_reactions)!=0]
CTL_genes_from_reactions=recon3D_gpr_rules_filtered[match(CTL_reactions, recon3D_gpr_rules_filtered[,1]),2]
CTL_genes_from_reactions = unique(unlist(strsplit(unlist(strsplit(gsub('[)]','', gsub('[(]','', na.omit(CTL_genes_from_reactions))), ' or ')), ' and ')))
CTL_genes_from_reactions = CTL_gene_PA_calls$final_gene_calls[match(CTL_genes_from_reactions, rownames(CTL_gene_PA_calls$final_gene_calls)),]
CTL_genes_from_reactions = rownames(CTL_genes_from_reactions)[CTL_genes_from_reactions[,'Final']!=0]
CTL_genes_from_reactions = recon3DModel_genes$ensembl[match(CTL_genes_from_reactions, recon3DModel_genes$entrez)]
# Universe of genes: all genes in recon3D model:
universe_genes = unique(recon3DModel_genes$ensembl)
# GO analysis:
CTL_GO_analysis_reactions = GO_analysis(CTL_genes_from_reactions, universe_genes)
CTL_GO_tables_reactions_sig = significant_GOs(CTL_GO_analysis_reactions)
View(CTL_GO_tables_reactions_sig$BP)
View(CTL_GO_tables_reactions_sig$CC)
View(CTL_GO_tables_reactions_sig$MF)
