#Directories:
base_dir = getwd() # When oppening with ptoject
Omics_data_Repos = '~/Documents/PhD/Omics_Data_Repos'

Metabolic_Models_Repo_general_code = paste(base_dir, 'general_code', sep='/')
Metabolic_Models_Repo_general_code_utils = paste(Metabolic_Models_Repo_general_code, 'utils', sep='/')
Metabolic_Models_Repo_general_code_R = paste(Metabolic_Models_Repo_general_code, 'R', sep='/')

Metabolic_Models_Repo_reconstruction_process_T_CELLS = paste(base_dir, 'pipeline_T_CELLS', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_Th1 = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS, 'Th1', sep='/')
Metabolic_Models_Repo_reconstruction_process_T_CELLS_Th1_PA_calls = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Th1, 'PA_calls', sep='/')

T_Cell_Data_Repo = paste(Omics_data_Repos, 'T_Cell_Data_Repo', sep='/')
T_Cell_Data_Repo_Data_Files = paste(T_Cell_Data_Repo, 'Data_Files', sep='/')
T_Cell_Data_Repo_Data_Files_Th1 = paste(T_Cell_Data_Repo_Data_Files, 'Th1', sep='/')
T_Cell_Data_Repo_Data_Files_Th1_RNA_seq = paste(T_Cell_Data_Repo_Data_Files_Th1, 'RNA_seq', sep='/')
T_Cell_Data_Repo_Data_Files_Th1_Microarray = paste(T_Cell_Data_Repo_Data_Files_Th1, 'Microarray', sep='/')
T_Cell_Data_Repo_Data_Files_Th1_Proteomics = paste(T_Cell_Data_Repo_Data_Files_Th1, 'Proteomics', sep='/')

#Files:
PA_calls_file = paste(Metabolic_Models_Repo_general_code_R, 'PA_calls.R', sep='/')
recon3D_consistent_genes_file = paste(Metabolic_Models_Repo_general_code_utils, 'entrez_genes/recon3D_consistent_genes.txt', sep='/')
recon3D_consistent_GPR_file = paste(Metabolic_Models_Repo_general_code_utils, 'GPRs/recon3D_consistent_GPR.txt', sep='/')
recon3DModel_gene_mapping_file = paste(Metabolic_Models_Repo_general_code_utils, 'recon3DModel_gene_mapping.csv', sep='/')

Th1_microarray_data_file = paste(T_Cell_Data_Repo_Data_Files_Th1_Microarray, 'Th1_microarray_data.csv', sep='/')
Th1_rnaseq_data_file = paste(T_Cell_Data_Repo_Data_Files_Th1_RNA_seq, 'Th1_rnaseq_data.csv', sep='/')
Th1_proteomics_data_file = paste(T_Cell_Data_Repo_Data_Files_Th1_Proteomics, 'Th1_proteomics_data.csv', sep='/')
Th1_microarray_metadata_file = paste(T_Cell_Data_Repo_Data_Files_Th1_Microarray, 'Th1_microarray_metadata.csv', sep='/')
Th1_rnaseq_metadata_file = paste(T_Cell_Data_Repo_Data_Files_Th1_RNA_seq, 'Th1_rnaseq_metadata.csv', sep='/')

#Source:
source(PA_calls_file)


#--- Read expression data for PA calls ---

recon3D_consistent_genes = as.character(read.table(recon3D_consistent_genes_file)[,1])

Th1_expression_data = read_omics_data(recon3D_consistent_genes, recon3DModel_gene_mapping_file,
                                      transcriptomics_data_files = list(RNAseq=Th1_rnaseq_data_file, Microarray=Th1_microarray_data_file),
                                      transcriptomics_metadata_files = list(RNAseq=Th1_rnaseq_metadata_file, Microarray=Th1_microarray_metadata_file),
                                      proteomics_data_files = list(ImmProt=Th1_proteomics_data_file),
                                      proteomics_metadata_files = NULL)

#--- Reaction Calls ---

Th1_PA_calls = PA_reaction_calls(Th1_expression_data, recon3DModel_gene_mapping_file, recon3D_consistent_GPR_file, recon3D_consistent_genes)

# Store data in files:
for(data_type in names(Th1_PA_calls$calls_per_sample$Transcriptomics)){
  file_name = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Th1_PA_calls, '/Th1_t_', data_type, '.csv', sep='')
  write.csv(Th1_PA_calls$calls_per_sample$Transcriptomics[[data_type]], file_name, row.names=T, col.names=T)
}
for(data_type in names(Th1_PA_calls$calls_per_sample$Proteomics)){
  file_name = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Th1_PA_calls, '/Th1_p_', data_type, '.csv', sep='')
  write.csv(Th1_PA_calls$calls_per_sample$Proteomics[[data_type]], file_name, row.names=T, col.names=T)
}
file_name = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Th1_PA_calls, 'Th1_gene_calls.csv', sep='/')
write.csv(Th1_PA_calls$gene_scores, file_name, row.names=T, col.names=F)
file_name = paste(Metabolic_Models_Repo_reconstruction_process_T_CELLS_Th1_PA_calls, 'Th1_reaction_calls.csv', sep='/')
write.csv(Th1_PA_calls$reaction_calls, file_name, row.names=T, col.names=T)



