#Directories:
base_dir = getwd() # When opening with project
Omics_data_Repos = '~/Documents/PhD/Omics_Data_Repos'

Metabolic_Models_Repo_general = paste(base_dir, 'GENERAL', sep='/')
Metabolic_Models_Repo_general_code = paste(Metabolic_Models_Repo_general, 'code', sep='/')
Metabolic_Models_Repo_general_code_R = paste(Metabolic_Models_Repo_general_code, 'R', sep='/')
Metabolic_Models_Repo_general_utility_data = paste(Metabolic_Models_Repo_general, 'utility_data', sep='/')

Metabolic_Models_Repo_model_reconstructions_T_CELLS = paste(base_dir, 'MODEL_RECONSTRUCTIONS', 'T_CELLS', sep='/')
Metabolic_Models_Repo_model_reconstructions_T_CELLS_Tfollicular = paste(Metabolic_Models_Repo_model_reconstructions_T_CELLS, 'subtypes_pipelines', 'Tfollicular', sep='/')
Metabolic_Models_Repo_model_reconstructions_T_CELLS_Tfollicular_PA_calls = paste(Metabolic_Models_Repo_model_reconstructions_T_CELLS_Tfollicular, 'PA_calls', sep='/')

T_Cell_Data_Repo = paste(Omics_data_Repos, 'T_Cell_Data_Repo', sep='/')
T_Cell_Data_Repo_Data_Files = paste(T_Cell_Data_Repo, 'Data_Files', sep='/')
T_Cell_Data_Repo_Data_Files_Tfollicular = paste(T_Cell_Data_Repo_Data_Files, 'Tfollicular', sep='/')
T_Cell_Data_Repo_Data_Files_Tfollicular_RNA_seq = paste(T_Cell_Data_Repo_Data_Files_Tfollicular, 'RNA_seq', sep='/')
T_Cell_Data_Repo_Data_Files_Tfollicular_Microarray = paste(T_Cell_Data_Repo_Data_Files_Tfollicular, 'Microarray', sep='/')

#Files:
PA_calls_file = paste(Metabolic_Models_Repo_general_code_R, 'PA_calls.R', sep='/')
HumanGEM_genes_file = paste(Metabolic_Models_Repo_general_utility_data, 'HumanGEM1.4.1_geneID_mapping.tsv', sep='/')
HumanGEM_GPR_file = paste(Metabolic_Models_Repo_general_utility_data, 'HumanGEM-1.4.1_GPRs.txt', sep='/')

Tfollicular_microarray_data_file = paste(T_Cell_Data_Repo_Data_Files_Tfollicular_Microarray, 'Tfollicular_microarray_data.csv', sep='/')
Tfollicular_rnaseq_data_file = paste(T_Cell_Data_Repo_Data_Files_Tfollicular_RNA_seq, 'Tfollicular_rnaseq_data.csv', sep='/')
Tfollicular_microarray_metadata_file = paste(T_Cell_Data_Repo_Data_Files_Tfollicular_Microarray, 'Tfollicular_microarray_metadata.csv', sep='/')
Tfollicular_rnaseq_metadata_file = paste(T_Cell_Data_Repo_Data_Files_Tfollicular_RNA_seq, 'Tfollicular_rnaseq_metadata.csv', sep='/')

#Source:
source(PA_calls_file)


#--- Read expression data for PA calls ---

HumanGEM_genes = as.character(unique(read.table(HumanGEM_genes_file, header=T, skip=1, sep='\t')$Gene_stable_ID))

Tfollicular_expression_data = read_omics_data(HumanGEM_genes,
                                       transcriptomics_data_files = list(RNAseq=Tfollicular_rnaseq_data_file, Microarray=Tfollicular_microarray_data_file),
                                       transcriptomics_metadata_files = list(RNAseq=Tfollicular_rnaseq_metadata_file, Microarray=Tfollicular_microarray_metadata_file),
                                       proteomics_data_files = NULL,
                                       proteomics_metadata_files = NULL)

#--- Reaction Calls ---

Tfollicular_PA_calls = PA_reaction_calls(Tfollicular_expression_data, HumanGEM_GPR_file, HumanGEM_genes)

# Store data in files:
for(data_type in names(Tfollicular_PA_calls$calls_per_sample$Transcriptomics)){
  file_name = paste(Metabolic_Models_Repo_model_reconstructions_T_CELLS_Tfollicular_PA_calls, '/Tfollicular_t_', data_type, '.csv', sep='')
  write.csv(Tfollicular_PA_calls$calls_per_sample$Transcriptomics[[data_type]], file_name, row.names=T, col.names=T)
}
file_name = paste(Metabolic_Models_Repo_model_reconstructions_T_CELLS_Tfollicular_PA_calls, 'Tfollicular_gene_calls.csv', sep='/')
write.csv(Tfollicular_PA_calls$gene_scores, file_name, row.names=T, col.names=F)
file_name = paste(Metabolic_Models_Repo_model_reconstructions_T_CELLS_Tfollicular_PA_calls, 'Tfollicular_reaction_calls.csv', sep='/')
write.csv(Tfollicular_PA_calls$reaction_calls, file_name, row.names=T, col.names=T)



