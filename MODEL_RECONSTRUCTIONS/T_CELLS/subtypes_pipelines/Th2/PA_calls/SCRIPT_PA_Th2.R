#Directories:
base_dir = getwd() # When opening with project
Omics_data_Repos = '~/Documents/PhD/Omics_Data_Repos'

Metabolic_Models_Repo_general = paste(base_dir, 'GENERAL', sep='/')
Metabolic_Models_Repo_general_code = paste(Metabolic_Models_Repo_general, 'code', sep='/')
Metabolic_Models_Repo_general_code_R = paste(Metabolic_Models_Repo_general_code, 'R', sep='/')
Metabolic_Models_Repo_general_utility_data = paste(Metabolic_Models_Repo_general, 'utility_data', sep='/')

Metabolic_Models_Repo_model_reconstructions_T_CELLS = paste(base_dir, 'MODEL_RECONSTRUCTIONS', 'T_CELLS', sep='/')
Metabolic_Models_Repo_model_reconstructions_T_CELLS_general_utility_data = paste(Metabolic_Models_Repo_model_reconstructions_T_CELLS, 'general/utility_data', sep='/')

Metabolic_Models_Repo_model_reconstructions_T_CELLS_Th2 = paste(Metabolic_Models_Repo_model_reconstructions_T_CELLS, 'subtypes_pipelines', 'Th2', sep='/')
Metabolic_Models_Repo_model_reconstructions_T_CELLS_Th2_PA_calls = paste(Metabolic_Models_Repo_model_reconstructions_T_CELLS_Th2, 'PA_calls', sep='/')

T_Cell_Data_Repo = paste(Omics_data_Repos, 'T_Cell_Data_Repo', sep='/')
T_Cell_Data_Repo_Data_Files = paste(T_Cell_Data_Repo, 'Data_Files', sep='/')
T_Cell_Data_Repo_Data_Files_Th2 = paste(T_Cell_Data_Repo_Data_Files, 'Th2', sep='/')
T_Cell_Data_Repo_Data_Files_Th2_RNA_seq = paste(T_Cell_Data_Repo_Data_Files_Th2, 'RNA_seq', sep='/')
T_Cell_Data_Repo_Data_Files_Th2_Microarray = paste(T_Cell_Data_Repo_Data_Files_Th2, 'Microarray', sep='/')
T_Cell_Data_Repo_Data_Files_Th2_Proteomics = paste(T_Cell_Data_Repo_Data_Files_Th2, 'Proteomics', sep='/')

#Files:
PA_calls_file = paste(Metabolic_Models_Repo_general_code_R, 'PA_calls.R', sep='/')
HumanGEM_genes_file = paste(Metabolic_Models_Repo_general_utility_data, 'HumanGEM-1.4.1_consistent_geneID_mapping.csv', sep='/')
HumanGEM_GPR_file = paste(Metabolic_Models_Repo_model_reconstructions_T_CELLS_general_utility_data, 'HumanGEM-1.4.1_forTcells_GPRs.txt', sep='/')

Th2_microarray_data_file = paste(T_Cell_Data_Repo_Data_Files_Th2_Microarray, 'Th2_microarray_data.csv', sep='/')
Th2_rnaseq_data_file = paste(T_Cell_Data_Repo_Data_Files_Th2_RNA_seq, 'Th2_rnaseq_data.csv', sep='/')
Th2_proteomics_data_file = paste(T_Cell_Data_Repo_Data_Files_Th2_Proteomics, 'Th2_proteomics_data.csv', sep='/')
Th2_microarray_metadata_file = paste(T_Cell_Data_Repo_Data_Files_Th2_Microarray, 'Th2_microarray_metadata.csv', sep='/')
Th2_rnaseq_metadata_file = paste(T_Cell_Data_Repo_Data_Files_Th2_RNA_seq, 'Th2_rnaseq_metadata.csv', sep='/')

#Source:
source(PA_calls_file)


#--- Read expression data for PA calls ---

HumanGEM_genes = as.character(unique(read.csv(HumanGEM_genes_file, header=T)$Gene_stable_ID))

Th2_expression_data = read_omics_data(HumanGEM_genes,
                                       transcriptomics_data_files = list(RNAseq=Th2_rnaseq_data_file, Microarray=Th2_microarray_data_file),
                                       transcriptomics_metadata_files = list(RNAseq=Th2_rnaseq_metadata_file, Microarray=Th2_microarray_metadata_file),
                                       proteomics_data_files = list(ImmProt=Th2_proteomics_data_file),
                                       proteomics_metadata_files = NULL)

#--- Reaction Calls ---

Th2_PA_calls = PA_reaction_calls(Th2_expression_data, HumanGEM_GPR_file, HumanGEM_genes)

# Store data in files:
for(data_type in names(Th2_PA_calls$calls_per_sample$Transcriptomics)){
  file_name = paste(Metabolic_Models_Repo_model_reconstructions_T_CELLS_Th2_PA_calls, '/Th2_t_', data_type, '.csv', sep='')
  write.csv(Th2_PA_calls$calls_per_sample$Transcriptomics[[data_type]], file_name, row.names=T, col.names=T)
}
for(data_type in names(Th2_PA_calls$calls_per_sample$Proteomics)){
  file_name = paste(Metabolic_Models_Repo_model_reconstructions_T_CELLS_Th2_PA_calls, '/Th2_p_', data_type, '.csv', sep='')
  write.csv(Th2_PA_calls$calls_per_sample$Proteomics[[data_type]], file_name, row.names=T, col.names=T)
}
file_name = paste(Metabolic_Models_Repo_model_reconstructions_T_CELLS_Th2_PA_calls, 'Th2_gene_calls.csv', sep='/')
write.csv(Th2_PA_calls$gene_scores, file_name, row.names=T, col.names=F)
file_name = paste(Metabolic_Models_Repo_model_reconstructions_T_CELLS_Th2_PA_calls, 'Th2_reaction_calls.csv', sep='/')
write.csv(Th2_PA_calls$reaction_calls, file_name, row.names=T, col.names=T)



