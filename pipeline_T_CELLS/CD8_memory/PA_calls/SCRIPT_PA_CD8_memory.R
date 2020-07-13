#Directories:
base_dir = getwd() # When oppening with ptoject
Omics_data_Repos = '~/Documents/PhD/Omics_Data_Repos'

Metabolic_Models_Repo_general_code = paste(base_dir, 'general_code', sep='/')
Metabolic_Models_Repo_general_code_utils = paste(Metabolic_Models_Repo_general_code, 'utils', sep='/')
Metabolic_Models_Repo_general_code_R = paste(Metabolic_Models_Repo_general_code, 'R', sep='/')

Metabolic_Models_Repo_reconstruction_process_T_CELLS = paste(base_dir, 'pipeline_T_CELLS', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_CD8_memory = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS, 'CD8_memory', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_CD8_memory_PA_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_CD8_memory, 'PA_calls', sep='/')

T_Cell_Data_Repo = paste(Omics_data_Repos, 'T_Cell_Data_Repo', sep='/')
T_Cell_Data_Repo_Data_Files = paste(T_Cell_Data_Repo, 'Data_Files', sep='/')
T_Cell_Data_Repo_Data_Files_CD8_memory = paste(T_Cell_Data_Repo_Data_Files, 'CD8_memory', sep='/')
T_Cell_Data_Repo_Data_Files_CD8_memory_RNA_seq = paste(T_Cell_Data_Repo_Data_Files_CD8_memory, 'RNA_seq', sep='/')
T_Cell_Data_Repo_Data_Files_CD8_memory_Microarray = paste(T_Cell_Data_Repo_Data_Files_CD8_memory, 'Microarray', sep='/')
T_Cell_Data_Repo_Data_Files_CD8_memory_Proteomics = paste(T_Cell_Data_Repo_Data_Files_CD8_memory, 'Proteomics', sep='/')

#Files:
PA_calls_file = paste(Metabolic_Models_Repo_general_code_R, 'PA_calls.R', sep='/')
recon3D_consistent_genes_file = paste(Metabolic_Models_Repo_general_code_utils, 'entrez_genes/recon3D_consistent_genes.txt', sep='/')
recon3D_consistent_GPR_file = paste(Metabolic_Models_Repo_general_code_utils, 'GPRs/recon3D_consistent_GPR.txt', sep='/')
recon3DModel_gene_mapping_file = paste(Metabolic_Models_Repo_general_code_utils, 'recon3DModel_gene_mapping.csv', sep='/')

CD8_memory_microarray_data_file = paste(T_Cell_Data_Repo_Data_Files_CD8_memory_Microarray, 'CD8_memory_microarray_data.csv', sep='/')
CD8_memory_rnaseq_data_file = paste(T_Cell_Data_Repo_Data_Files_CD8_memory_RNA_seq, 'CD8_memory_rnaseq_data.csv', sep='/')
CD8_memory_proteomics_data_file = paste(T_Cell_Data_Repo_Data_Files_CD8_memory_Proteomics, 'CD8_memory_proteomics_data.csv', sep='/')
CD8_memory_microarray_metadata_file = paste(T_Cell_Data_Repo_Data_Files_CD8_memory_Microarray, 'CD8_memory_microarray_metadata.csv', sep='/')
CD8_memory_rnaseq_metadata_file = paste(T_Cell_Data_Repo_Data_Files_CD8_memory_RNA_seq, 'CD8_memory_rnaseq_metadata.csv', sep='/')

#Source:
source(PA_calls_file)


#--- Read expression data for PA calls ---

recon3D_consistent_genes = as.character(read.table(recon3D_consistent_genes_file)[,1])

CD8_memory_expression_data = read_omics_data(recon3D_consistent_genes, recon3DModel_gene_mapping_file,
                                             transcriptomics_data_files = list(RNAseq=CD8_memory_rnaseq_data_file, Microarray=CD8_memory_microarray_data_file),
                                             transcriptomics_metadata_files = list(RNAseq=CD8_memory_rnaseq_metadata_file, Microarray=CD8_memory_microarray_metadata_file),
                                             proteomics_data_files = list(ImmProt=CD8_memory_proteomics_data_file),
                                             proteomics_metadata_files = NULL)

#--- Reaction Calls ---

CD8_memory_PA_calls = PA_reaction_calls(CD8_memory_expression_data, recon3DModel_gene_mapping_file, recon3D_consistent_GPR_file, recon3D_consistent_genes)

# Store data in files:
for(data_type in names(CD8_memory_PA_calls$calls_per_sample$Transcriptomics)){
  file_name = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_CD8_memory_PA_calls, '/CD8_memory_t_', data_type, '.csv', sep='')
  write.csv(CD8_memory_PA_calls$calls_per_sample$Transcriptomics[[data_type]], file_name, row.names=T, col.names=T)
}
for(data_type in names(CD8_memory_PA_calls$calls_per_sample$Proteomics)){
  file_name = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_CD8_memory_PA_calls, '/CD8_memory_p_', data_type, '.csv', sep='')
  write.csv(CD8_memory_PA_calls$calls_per_sample$Proteomics[[data_type]], file_name, row.names=T, col.names=T)
}
file_name = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_CD8_memory_PA_calls, 'CD8_memory_gene_calls.csv', sep='/')
write.csv(CD8_memory_PA_calls$gene_scores, file_name, row.names=T, col.names=F)
file_name = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_CD8_memory_PA_calls, 'CD8_memory_reaction_calls.csv', sep='/')
write.csv(CD8_memory_PA_calls$reaction_calls, file_name, row.names=T, col.names=T)



