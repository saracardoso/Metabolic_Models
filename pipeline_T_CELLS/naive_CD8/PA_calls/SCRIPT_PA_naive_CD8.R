#Directories:
base_dir = getwd() # When oppening with ptoject
Omics_data_Repos = '~/Documents/PhD/Omics_Data_Repos'

Metabolic_Models_Repo_general_code = paste(base_dir, 'general_code', sep='/')
Metabolic_Models_Repo_general_code_utils = paste(Metabolic_Models_Repo_general_code, 'utils', sep='/')
Metabolic_Models_Repo_general_code_R = paste(Metabolic_Models_Repo_general_code, 'R', sep='/')

Metabolic_Models_Repo_reconstruction_process_T_CELLS = paste(base_dir, 'pipeline_T_CELLS', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD8 = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS, 'naive_CD8', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD8_PA_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD8, 'PA_calls', sep='/')

T_Cell_Data_Repo = paste(Omics_data_Repos, 'T_Cell_Data_Repo', sep='/')
T_Cell_Data_Repo_Data_Files = paste(T_Cell_Data_Repo, 'Data_Files', sep='/')
T_Cell_Data_Repo_Data_Files_naive_CD8 = paste(T_Cell_Data_Repo_Data_Files, 'naive_CD8', sep='/')
T_Cell_Data_Repo_Data_Files_naive_CD8_RNA_seq = paste(T_Cell_Data_Repo_Data_Files_naive_CD8, 'RNA_seq', sep='/')
T_Cell_Data_Repo_Data_Files_naive_CD8_Microarray = paste(T_Cell_Data_Repo_Data_Files_naive_CD8, 'Microarray', sep='/')
T_Cell_Data_Repo_Data_Files_naive_CD8_Proteomics = paste(T_Cell_Data_Repo_Data_Files_naive_CD8, 'Proteomics', sep='/')

#Files:
PA_calls_file = paste(Metabolic_Models_Repo_general_code_R, 'PA_calls.R', sep='/')
recon3D_consistent_genes_file = paste(Metabolic_Models_Repo_general_code_utils, 'entrez_genes/recon3D_consistent_genes.txt', sep='/')
recon3D_consistent_GPR_file = paste(Metabolic_Models_Repo_general_code_utils, 'GPRs/recon3D_consistent_GPR.txt', sep='/')
recon3DModel_gene_mapping_file = paste(Metabolic_Models_Repo_general_code_utils, 'recon3DModel_gene_mapping.csv', sep='/')

naive_CD8_microarray_data_file = paste(T_Cell_Data_Repo_Data_Files_naive_CD8_Microarray, 'naive_CD8_microarray_data.csv', sep='/')
naive_CD8_rnaseq_data_file = paste(T_Cell_Data_Repo_Data_Files_naive_CD8_RNA_seq, 'naive_CD8_rnaseq_data.csv', sep='/')
naive_CD8_proteomics_data_file = paste(T_Cell_Data_Repo_Data_Files_naive_CD8_Proteomics, 'naive_CD8_proteomics_data.csv', sep='/')
naive_CD8_microarray_metadata_file = paste(T_Cell_Data_Repo_Data_Files_naive_CD8_Microarray, 'naive_CD8_microarray_metadata.csv', sep='/')
naive_CD8_rnaseq_metadata_file = paste(T_Cell_Data_Repo_Data_Files_naive_CD8_RNA_seq, 'naive_CD8_rnaseq_metadata.csv', sep='/')

#Source:
source(PA_calls_file)


#--- Read expression data for PA calls ---

recon3D_consistent_genes = as.character(read.table(recon3D_consistent_genes_file)[,1])

naive_CD8_expression_data = read_omics_data(recon3D_consistent_genes, recon3DModel_gene_mapping_file,
                                       transcriptomics_data_files = list(RNAseq=naive_CD8_rnaseq_data_file, Microarray=naive_CD8_microarray_data_file),
                                       transcriptomics_metadata_files = list(RNAseq=naive_CD8_rnaseq_metadata_file, Microarray=naive_CD8_microarray_metadata_file),
                                       proteomics_data_files = list(ImmProt=naive_CD8_proteomics_data_file),
                                       proteomics_metadata_files = NULL)

#--- Reaction Calls ---

naive_CD8_PA_calls = PA_reaction_calls(naive_CD8_expression_data, recon3DModel_gene_mapping_file, recon3D_consistent_GPR_file, recon3D_consistent_genes)

# Store data in files:
for(data_type in names(naive_CD8_PA_calls$calls_per_sample$Transcriptomics)){
  file_name = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD8_PA_calls, '/naive_CD8_t_', data_type, '.csv', sep='')
  write.csv(naive_CD8_PA_calls$calls_per_sample$Transcriptomics[[data_type]], file_name, row.names=T, col.names=T)
}
for(data_type in names(naive_CD8_PA_calls$calls_per_sample$Proteomics)){
  file_name = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD8_PA_calls, '/naive_CD8_p_', data_type, '.csv', sep='')
  write.csv(naive_CD8_PA_calls$calls_per_sample$Proteomics[[data_type]], file_name, row.names=T, col.names=T)
}
file_name = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD8_PA_calls, 'naive_CD8_gene_calls.csv', sep='/')
write.csv(naive_CD8_PA_calls$gene_scores, file_name, row.names=T, col.names=F)
file_name = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_naive_CD8_PA_calls, 'naive_CD8_reaction_calls.csv', sep='/')
write.csv(naive_CD8_PA_calls$reaction_calls, file_name, row.names=T, col.names=T)



