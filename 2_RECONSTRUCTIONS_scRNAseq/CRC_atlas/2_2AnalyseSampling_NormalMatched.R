code_dir = './GENERAL/code/R'
for(file in list.files(code_dir, full.names=TRUE)) source(file)

HumanGEM_rxns_subsystems = jsonlite::read_json('./GENERAL/utility_data/reactions_subsystems_mapping.json',
                                               simplifyVector = TRUE)
HumanGEM_subsystems_rxns = jsonlite::read_json('./GENERAL/utility_data/subsystems_reactions_mapping.json',
                                               simplifyVector = TRUE)

ct_colors = c('#b30000', '#999999', '#cc79a7', '#b3c9e6', '#4e84c4', '#c3d7a4', '#52854c', '#fdb981',
              '#d16103', '#293352')
names(ct_colors) = c('Cytotoxic CD8 Tcells', 'Follicular CD4 Tcells', 'IL17+ CD4 Tcells',
                     'Memory CD4 Tcells', 'Memory CD8 Tcells', 'Naive CD4 Tcells', 'Naive CD8 Tcells',
                     'Proliferative CD4 Tcells', 'Proliferative CD8 Tcells',  'Regulatory CD4 Tcells')

state_colors = c('#52854c', '#b30000', '#ff6666', '#cc79a7')
names(state_colors) = c('Normal Matched', 'Tumour', 'Tumour Core', 'Tumour Border')

individual_colors = c('#b30000', '#999999', '#cc79a7', '#b3c9e6', '#4e84c4', '#c3d7a4', '#52854c',
                      '#fdb981', '#d16103', '#293352', '#ffdb6d', '#c4961a', '#003366')
names(individual_colors) = c('31', '32', '33', '35', 'KUL01', 'KUL19', 'KUL21', 'SMC01', 'SMC04', 'SMC06',
                             'SMC07', 'SMC08', 'SMC10')



# ---------------
# --- CONTROL ---
# ---------------



# -----
# --- Create json file with names of models that 'passed' the sampling analysis [Control]
# --- [all samplings will be stored in this file]
# -----

CRCatlas_sampling = list()
CRCatlas_sampling$NormalMatched = list()

# Patient 31
CRCatlas_sampling$NormalMatched[['31']] = list()
CRCatlas_sampling$NormalMatched[['31']]$control = list()
CRCatlas_sampling$NormalMatched[['31']]$control$scrEXT001 = c('IL17+ CD4 Tcells', 'Memory CD4 Tcells',
                                                              'Proliferative CD4 Tcells',
                                                              'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['31']]$control$scrEXT002 = c('Cytotoxic CD8 Tcells',
                                                              'Memory CD4 Tcells',
                                                              'Naive CD4 Tcells',
                                                              'Proliferative CD8 Tcells',
                                                              'Regulatory CD4 Tcells',
                                                              'Follicular CD4 Tcells')
CRCatlas_sampling$NormalMatched[['31']]$control$scrEXT003 = c('Memory CD4 Tcells', 'Naive CD4 Tcells')

# Patient 32
CRCatlas_sampling$NormalMatched[['32']] = list()
CRCatlas_sampling$NormalMatched[['32']]$control = list()
CRCatlas_sampling$NormalMatched[['32']]$control$scrEXT009 = c('Cytotoxic CD8 Tcells',
                                                              'Follicular CD4 Tcells',
                                                              'IL17+ CD4 Tcells', 'Memory CD4 Tcells',
                                                              'Memory CD8 Tcells', 'Naive CD4 Tcells',
                                                              'Naive CD8 Tcells',
                                                              'Proliferative CD4 Tcells',
                                                              'Proliferative CD8 Tcells',
                                                              'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['32']]$control$scrEXT010 = c('Cytotoxic CD8 Tcells',
                                                              'Follicular CD4 Tcells',
                                                              'Memory CD4 Tcells',
                                                              'Memory CD8 Tcells', 'Naive CD4 Tcells',
                                                              'Naive CD8 Tcells',
                                                              'Proliferative CD4 Tcells',
                                                              'Proliferative CD8 Tcells',
                                                              'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['32']]$control$scrEXT011 = c('Cytotoxic CD8 Tcells',
                                                              'Memory CD8 Tcells', 'Naive CD4 Tcells',
                                                              'Regulatory CD4 Tcells')

# Patient 33
CRCatlas_sampling$NormalMatched[['33']] = list()
CRCatlas_sampling$NormalMatched[['33']]$control = list()
CRCatlas_sampling$NormalMatched[['33']]$control$scrEXT012 = c('Memory CD4 Tcells', 'Memory CD8 Tcells',
                                                              'Proliferative CD4 Tcells',
                                                              'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['33']]$control$scrEXT013 = c('Memory CD4 Tcells',
                                                              'Naive CD4 Tcells',
                                                              'Proliferative CD4 Tcells',
                                                              'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['33']]$control$scrEXT014 = c('IL17+ CD4 Tcells',
                                                              'Memory CD4 Tcells',
                                                              'Memory CD8 Tcells', 'Naive CD4 Tcells')

# Patient 35
CRCatlas_sampling$NormalMatched[['35']] = list()
CRCatlas_sampling$NormalMatched[['35']]$control = list()
CRCatlas_sampling$NormalMatched[['35']]$control$scrEXT018 = c('Follicular CD4 Tcells',
                                                              'IL17+ CD4 Tcells', 'Memory CD4 Tcells',
                                                              'Memory CD8 Tcells',
                                                              'Naive CD4 Tcells',
                                                              'Proliferative CD4 Tcells',
                                                              'Proliferative CD8 Tcells')
CRCatlas_sampling$NormalMatched[['35']]$control$scrEXT019 = c('Follicular CD4 Tcells',
                                                              'Memory CD4 Tcells',
                                                              'Memory CD8 Tcells',
                                                              'Naive CD4 Tcells',
                                                              'Proliferative CD4 Tcells',
                                                              'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['35']]$control$scrEXT020 = c('Memory CD8 Tcells', 'Naive CD4 Tcells',
                                                              'Regulatory CD4 Tcells')

# Patient KUL01
CRCatlas_sampling$NormalMatched[['KUL01']] = list()
CRCatlas_sampling$NormalMatched[['KUL01']]$control = list()
CRCatlas_sampling$NormalMatched[['KUL01']]$control[['KUL01-B']] = c('Cytotoxic CD8 Tcells',
                                                                    'Follicular CD4 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['KUL01']]$control[['KUL01-N']] = c('Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells')
CRCatlas_sampling$NormalMatched[['KUL01']]$control[['KUL01-T']] = c('IL17+ CD4 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells')

# Patient KUL19
CRCatlas_sampling$NormalMatched[['KUL19']] = list()
CRCatlas_sampling$NormalMatched[['KUL19']]$control = list()
CRCatlas_sampling$NormalMatched[['KUL19']]$control[['KUL19-B']] = c('Cytotoxic CD8 Tcells',
                                                                    'Follicular CD4 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Naive CD8 Tcells',
                                                                    'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['KUL19']]$control[['KUL19-N']] = c('Cytotoxic CD8 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['KUL19']]$control[['KUL19-T']] = c('Cytotoxic CD8 Tcells',
                                                                    'Follicular CD4 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Naive CD8 Tcells',
                                                                    'Regulatory CD4 Tcells')

# Patient KUL21
CRCatlas_sampling$NormalMatched[['KUL21']] = list()
CRCatlas_sampling$NormalMatched[['KUL21']]$control = list()
CRCatlas_sampling$NormalMatched[['KUL21']]$control[['KUL21-B']] = c('Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['KUL21']]$control[['KUL21-N']] = c('IL17+ CD4 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['KUL21']]$control[['KUL21-T']] = c('Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')

# Patient SMC01
CRCatlas_sampling$NormalMatched[['SMC01']] = list()
CRCatlas_sampling$NormalMatched[['SMC01']]$control = list()
CRCatlas_sampling$NormalMatched[['SMC01']]$control[['SMC01-N']] = c('Cytotoxic CD8 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['SMC01']]$control[['SMC01-T']] = c('Cytotoxic CD8 Tcells',
                                                                    'Follicular CD4 Tcells',
                                                                    'IL17+ CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Proliferative CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')

# Patient SMC04
CRCatlas_sampling$NormalMatched[['SMC04']] = list()
CRCatlas_sampling$NormalMatched[['SMC04']]$control = list()
CRCatlas_sampling$NormalMatched[['SMC04']]$control[['SMC04-N']] = c('Cytotoxic CD8 Tcells',
                                                                    'IL17+ CD4 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Naive CD8 Tcells')
CRCatlas_sampling$NormalMatched[['SMC04']]$control[['SMC04-T']] = c('Cytotoxic CD8 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Naive CD8 Tcells',
                                                                    'Regulatory CD4 Tcells')

# Patient SMC06
CRCatlas_sampling$NormalMatched[['SMC06']] = list()
CRCatlas_sampling$NormalMatched[['SMC06']]$control = list()
CRCatlas_sampling$NormalMatched[['SMC06']]$control[['SMC06-N']] = c('Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Naive CD8 Tcells',
                                                                    'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['SMC06']]$control[['SMC06-T']] = c('Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')

# Patient SMC07
CRCatlas_sampling$NormalMatched[['SMC07']] = list()
CRCatlas_sampling$NormalMatched[['SMC07']]$control = list()
CRCatlas_sampling$NormalMatched[['SMC07']]$control[['SMC07-N']] = c('Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['SMC07']]$control[['SMC07-T']] = c('Cytotoxic CD8 Tcells',
                                                                    'Follicular CD4 Tcells',
                                                                    'IL17+ CD4 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Naive CD8 Tcells',
                                                                    'Proliferative CD4 Tcells',
                                                                    'Proliferative CD8 Tcells',
                                                                    'Regulatory CD4 Tcells')

# Patient SMC08
CRCatlas_sampling$NormalMatched[['SMC08']] = list()
CRCatlas_sampling$NormalMatched[['SMC08']]$control = list()
CRCatlas_sampling$NormalMatched[['SMC08']]$control[['SMC08-N']] = c('Cytotoxic CD8 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['SMC08']]$control[['SMC08-T']] = c('Cytotoxic CD8 Tcells',
                                                                    'Follicular CD4 Tcells',
                                                                    'IL17+ CD4 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Proliferative CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')

# Patient SMC10
CRCatlas_sampling$NormalMatched[['SMC10']] = list()
CRCatlas_sampling$NormalMatched[['SMC10']]$control = list()
CRCatlas_sampling$NormalMatched[['SMC10']]$control[['SMC10-N']] = c('Cytotoxic CD8 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['SMC10']]$control[['SMC10-T']] = c('Cytotoxic CD8 Tcells',
                                                                    'Follicular CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Proliferative CD4 Tcells',
                                                                    'Proliferative CD8 Tcells',
                                                                    'Regulatory CD4 Tcells')

# Save the list into a json file:
jsonlite::write_json(CRCatlas_sampling, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/CRCatlas_sampling.json')





# -----
# --- Read sampling data and metadata
# -----

# Read metadata of samples:
metadata = read.csv('./0Data/scRNAseq/CRC_atlas/expression_data/metadata.csv', row.names=1)

# Read sampling results and organising sampling metadata into one big dataframe:
samples = c()
individuals = c()
states = c()
cell_types = c()
media = c()
row_names = c()
CRCatlas_samplingNMControl_matrix = c()
CRCatlas_sampling_dataframe = data.table::melt(CRCatlas_sampling)
for(idx in 1:dim(CRCatlas_sampling_dataframe)[1]){
  individual = CRCatlas_sampling_dataframe[idx, 4]
  sample = CRCatlas_sampling_dataframe[idx, 2]
  cell_type = CRCatlas_sampling_dataframe[idx, 1]
  
  message(paste(CRCatlas_sampling_dataframe[idx,], collapse='_'), 
          ' ', idx, '/', dim(CRCatlas_sampling_dataframe)[1])
  
  # Sampling data:
  for(medium_type in c('Blood_SMDB')){
    # Add to list:
    message('| Reading data from medium ', medium_type, '...')
    mtx = Matrix::Matrix(as.matrix(read.csv(paste('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched',
                                                  individual, '1_sampling/control', sample, gsub(' ', '_', cell_type),
                                                  paste(medium_type, '.csv', sep=''), sep='/'), row.names=1)), sparse=T)
    # Change rownames:
    message('- Changing row names')
    rownms = paste(sample, gsub(' ', '_', cell_type), medium_type, rownames(mtx), sep='_')
    rownames(mtx) = rownms
    
    # Randomly choose only 500 samplings:
    message('- Choosing randomly 500 samplings...')
    samplings = sample(rownms, 500)
    mtx = mtx[samplings,]
    invisible(gc())
    
    # Add to matrix:
    message('- Adding new data to the general matrix...')
    CRCatlas_samplingNMControl_matrix = rbind(CRCatlas_samplingNMControl_matrix, mtx)
    
    # Save data until now:
    message('- Saving the matrix we have until now...')
    Matrix::writeMM(CRCatlas_samplingNMControl_matrix,
                    './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/sampling_control_data.mtx')
    
    # Samples metadata: [sample, individual, state, cell_type]
    message('- Collecting metadata...')
    n_samplings = dim(mtx)[1]
    samples = c(samples, rep(sample, n_samplings))
    individuals = c(individuals, rep(individual, n_samplings))
    states = c(states, rep(metadata[sample, 'Sample.Source'], n_samplings))
    cell_types = c(cell_types, rep(cell_type, n_samplings))
    media = c(media, rep(medium_type, n_samplings))
    row_names = c(row_names, samplings)
    
    # Deleting temporary mtx:
    message('- Deleting temporary mtx...')
    remove(mtx)
    invisible(gc())
  }
  message('')
}
CRCatlas_samplingNMControl_meta = data.frame(sample=samples, individual=individuals, state=states,
                                             cell_type=cell_types, medium=media, row.names=row_names)
jsonlite::write_json(dimnames(CRCatlas_samplingNMControl_matrix),
                     './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/sampling_control_data_dimnames.json')

# Save metadata:
write.csv(CRCatlas_samplingNMControl_meta,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/sampling_control_metadata.csv')





# -----
# --- Compare presence of reactions
# -----

# Read reactions presence file and metadata
reactions_presence = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/reaction_presence.csv',
                              row.names=1, check.names=FALSE)
#split_colnames = strsplit(colnames(reactions_presence), '_')
#for(idx in 1:length(split_colnames)){
#  colnames(reactions_presence)[idx] = gsub('[ ]', '_', paste(c(split_colnames[[idx]][2:3], 'Blood_SMDB'), collapse='_'))
#}
#write.csv(reactions_presence, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/reaction_presence.csv')
split_colnames = strsplit(colnames(reactions_presence), '_')
samples = c()
cell_types = c()
media = c()
for(row_split in split_colnames){
  samples = c(samples, row_split[1])
  cell_types = c(cell_types, paste(row_split[2:(length(split_colnames[[1]])-2)], collapse=' '))
  media = c(media, paste(tail(split_colnames[[1]], n=2), collapse='_'))
}
reactionsPresence_metadata = data.frame(sample=samples, cell_type=cell_types, medium=media,
                                        row.names=colnames(reactions_presence))
more_meta = read.csv('./0Data/scRNAseq/CRC_atlas/expression_data/metadata.csv', row.names=1)
states = c()
for(samp in reactionsPresence_metadata$sample){
  states = c(states, more_meta[samp, 'Sample.Source'])
}
reactionsPresence_metadata$state = states


nRxns_perModel = colSums(reactions_presence)
# Average number of reactions for each cell-type
for(ct in unique(reactionsPresence_metadata$cell_type)){
  message(ct)
  print(mean(nRxns_perModel[reactionsPresence_metadata$cell_type==ct]))
}

# Average number of reactions for each cell-type in each sample state
for(ct in unique(reactionsPresence_metadata$cell_type)){
  for(state in unique(reactionsPresence_metadata$state)){
    message(ct, ' - ', state)
    print(mean(nRxns_perModel[reactionsPresence_metadata$cell_type==ct & 
                                reactionsPresence_metadata$state==state]))
  }
}

# % of reaction overlapping for models in the same cell-type
message('Reactions overlap for each cell-type')
for(ct in unique(reactionsPresence_metadata$cell_type)){
  n_rxns = dim(reactions_presence)[1]
  n_models_ct = sum(reactionsPresence_metadata$cell_type==ct)
  overlapRXNS = rowSums(reactions_presence[,reactionsPresence_metadata$cell_type==ct])
  message(ct, ': ', sum(overlapRXNS == n_models_ct | overlapRXNS == 0) / n_rxns * 100,
          '(', n_models_ct, ')')
}
# What happens to these values if we allow for an interval (only 90% overlapp to consider a reaction
# present in 'all' models)??
message('Reactions overlap (50%) for each cell-type')
for(ct in unique(reactionsPresence_metadata$cell_type)){
  n_rxns = dim(reactions_presence)[1]
  n_models_ct = sum(reactionsPresence_metadata$cell_type==ct)
  overlapRXNS = rowSums(reactions_presence[,reactionsPresence_metadata$cell_type==ct])
  message(ct, ': ', sum(overlapRXNS >= .5*n_models_ct | overlapRXNS < .1*n_models_ct) / n_rxns * 100,
          '(', n_models_ct, ')')
}

# Does this mean that models of the same cell-type won't be able to group together based on reaction
# presence?? Checked next:

# PCA of reaction presence -> colour by: cell-type, state, (media)
pca.results = irlba::irlba(A=t(reactions_presence), nv=50)
rxn_embeddings = pca.results$u %*% diag(pca.results$d)

pca_plot_df = cbind(rxn_embeddings[,1:2], reactionsPresence_metadata)
colnames(pca_plot_df)[1:2] = c('PCA1', 'PCA2')
plt = ggplot2::ggplot(pca_plot_df, ggplot2::aes(PCA1, PCA2))
plt + ggplot2::geom_point(ggplot2::aes_string(colour='cell_type'))
plt + ggplot2::geom_point(ggplot2::aes_string(colour='medium'))
plt + ggplot2::geom_point(ggplot2::aes_string(colour='state')) #####
plt + ggplot2::geom_point(ggplot2::aes_string(colour='sample'))

# Could we see a separation in the PCA if we clustered the samplings?



# -----
# --- Analyse fluxes of reactions
# -----

# 1. Read data:
sampling_control_data = Matrix::readMM('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/sampling_control_data.mtx')
invisible(gc())
dimnames_data = jsonlite::read_json('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/sampling_control_data_dimnames.json',
                                       simplifyVector = TRUE)
dimnames(sampling_control_data) = dimnames_data
sampling_control_meta = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/sampling_control_metadata.csv',
                                 row.names=1, stringsAsFactors=TRUE)

# 2. Sigmoid transform the data:
sampling_control_data = sigmoid_norm(sampling_control_data)
invisible(gc())

# 3. Create SeuratObject, it eases performing certain analyses:
seurat_control = SeuratObject::CreateSeuratObject(Matrix::t(sampling_control_data), project='samplings',
                                                assay='RNA', meta.data = sampling_control_meta)
remove(sampling_control_data, sampling_control_meta)
invisible(gc())
# 3.1. Save the seurat dataset until now:
SeuratDisk::SaveH5Seurat(seurat_control, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/control.h5Seurat')

# 4. Calculate top 2000 variable reactions:
seurat_control = Seurat::FindVariableFeatures(seurat_control)
# 4.1. Plot variable reactions:
varPlot1 = Seurat::VariableFeaturePlot(seurat_control)
varPlot2 = Seurat::LabelPoints(varPlot1, points = head(Seurat::VariableFeatures(seurat_control), 10),
                               repel = TRUE)
varPlot2 + ggplot2::ylim(0, 2e32) + ggplot2::xlab('Average Flux')
# 4.2. Top 20 reactions:
head(Seurat::VariableFeatures(seurat_control), 20)
# MAR00611, MAR01656, MAR09254, MAR05111, MAR07753, MAR09062, MAR00444, MAR00639, MAR08838, MAR07726,
# MAR01659, MAR08837, MAR00583, MAR09873, MAR05525, MAR04567, MAR05527, MAR07190, MAR05182, MAR08828

# 5. Put our values in the scale.data slot, so that we can run the next analyses
seurat_control = Seurat::SetAssayData(seurat_control, assay='RNA', slot='scale.data',
                                      new.data=as.matrix(seurat_control@assays$RNA@counts))
seurat_control@assays$RNA@scale.data = seurat_control@assays$RNA@scale.data[Seurat::VariableFeatures(seurat_control),]
invisible(gc())
# 5.1. Save the seurat dataset until now:
SeuratDisk::SaveH5Seurat(seurat_control, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/control.h5Seurat',
                         overwrite=TRUE)
invisible(gc())

# 6. PCA on all models together:
seurat_control = Seurat::RunPCA(seurat_control, features = Seurat::VariableFeatures(seurat_control))
invisible(gc())
# 6.2. PCA plots coloured by cell-types, state amd individuals
Seurat::DimPlot(seurat_control, reduction='pca', group.by='cell_type', cols=ct_colors)
Seurat::DimPlot(seurat_control, reduction='pca', group.by='state', cols=state_colors)
Seurat::DimPlot(seurat_control, reduction='pca', group.by='individual', cols=individual_colors)
# 6.3. Check UMAP
# 6.3.1. Calculate the elbow:
pct = seurat_control[["pca"]]@stdev /sum(seurat_control[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
elbow = min(co1, co2) # elbow is 10
plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
  ggplot2::geom_text() + 
  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
# 6.3.2. Run UMAP
seurat_control = Seurat::RunUMAP(seurat_control, dims=1:elbow, reduction='pca')
# 6.3.3. Visualize UMAP
Seurat::DimPlot(seurat_control, reduction="umap", group.by='cell_type', label=FALSE, cols=ct_colors)
Seurat::DimPlot(seurat_control, reduction="umap", group.by='state', label=FALSE, cols=state_colors)
Seurat::DimPlot(seurat_control, reduction="umap", group.by='individual', label=FALSE,
                cols=individual_colors)
# 6.4. Save Seurat object:
SeuratDisk::SaveH5Seurat(seurat_control, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/control.h5Seurat',
                         overwrite=TRUE)
invisible(gc())
# 6.5. Save dataframe with PCA and UMAP information
infoDF = data.frame(PC1=seurat_control@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control@reductions$pca@cell.embeddings[,'PC_2'],
                    UMAP1=seurat_control@reductions$umap@cell.embeddings[,'UMAP_1'],
                    UMAP2=seurat_control@reductions$umap@cell.embeddings[,'UMAP_2'],
                    cell_type=seurat_control$cell_type, state=seurat_control$state,
                    individual=seurat_control$individual)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/all_samplings_reductions.csv')


# 7. PCA for each individual separately
# 7.1. Individual 31
seurat_control_31 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$individual=='31'])
# 7.1.1. Run PCA
seurat_control_31 = Seurat::RunPCA(seurat_control_31, features = Seurat::VariableFeatures(seurat_control_31))
invisible(gc())
# 7.1.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_31, reduction='pca', group.by='cell_type', cols=ct_colors)
Seurat::DimPlot(seurat_control_31, reduction='pca', group.by='state', cols=state_colors)
# 7.1.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_31@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_31@reductions$pca@cell.embeddings[,'PC_2'],
                    cell_type=seurat_control_31$cell_type, state=seurat_control_31$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/31_samplings_reductions.csv')

# 7.2. Individual 32
seurat_control_32 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$individual=='32'])
invisible(gc())
# 7.2.1. Run PCA
seurat_control_32 = Seurat::RunPCA(seurat_control_32, features = Seurat::VariableFeatures(seurat_control_32))
invisible(gc())
# 7.2.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_32, reduction='pca', group.by='cell_type', cols=ct_colors) /
Seurat::DimPlot(seurat_control_32, reduction='pca', group.by='state', cols=state_colors)
# 7.2.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_32@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_32@reductions$pca@cell.embeddings[,'PC_2'],
                    cell_type=seurat_control_32$cell_type, state=seurat_control_32$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/32_samplings_reductions.csv')

# 7.3. Individual 33
seurat_control_33 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$individual=='33'])
invisible(gc())
# 7.3.1. Run PCA
seurat_control_33 = Seurat::RunPCA(seurat_control_33, features = Seurat::VariableFeatures(seurat_control_33))
invisible(gc())
# 7.3.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_33, reduction='pca', group.by='cell_type', cols=ct_colors) /
  Seurat::DimPlot(seurat_control_33, reduction='pca', group.by='state', cols=state_colors)
# 7.3.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_33@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_33@reductions$pca@cell.embeddings[,'PC_2'],
                    cell_type=seurat_control_33$cell_type, state=seurat_control_33$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/33_samplings_reductions.csv')

# 7.4. Individual 35
seurat_control_35 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$individual=='35'])
invisible(gc())
# 7.4.1. Run PCA
seurat_control_35 = Seurat::RunPCA(seurat_control_35, features = Seurat::VariableFeatures(seurat_control_35))
invisible(gc())
# 7.4.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_35, reduction='pca', group.by='cell_type', cols=ct_colors) /
  Seurat::DimPlot(seurat_control_35, reduction='pca', group.by='state', cols=state_colors)
# 7.4.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_35@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_35@reductions$pca@cell.embeddings[,'PC_2'],
                    cell_type=seurat_control_35$cell_type, state=seurat_control_35$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/35_samplings_reductions.csv')

# 7.5. Individual KUL01
seurat_control_KUL01 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$individual=='KUL01'])
invisible(gc())
# 7.5.1. Run PCA
seurat_control_KUL01 = Seurat::RunPCA(seurat_control_KUL01, features = Seurat::VariableFeatures(seurat_control_KUL01))
invisible(gc())
# 7.5.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_KUL01, reduction='pca', group.by='cell_type', cols=ct_colors) /
  Seurat::DimPlot(seurat_control_KUL01, reduction='pca', group.by='state', cols=state_colors)
# 7.5.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_KUL01@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_KUL01@reductions$pca@cell.embeddings[,'PC_2'],
                    cell_type=seurat_control_KUL01$cell_type, state=seurat_control_KUL01$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/KUL01_samplings_reductions.csv')

# 7.6. Individual KUL19
seurat_control_KUL19 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$individual=='KUL19'])
invisible(gc())
# 7.6.1. Run PCA
seurat_control_KUL19 = Seurat::RunPCA(seurat_control_KUL19, features = Seurat::VariableFeatures(seurat_control_KUL19))
invisible(gc())
# 7.6.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_KUL19, reduction='pca', group.by='cell_type', cols=ct_colors) /
  Seurat::DimPlot(seurat_control_KUL19, reduction='pca', group.by='state', cols=state_colors)
# 7.6.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_KUL19@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_KUL19@reductions$pca@cell.embeddings[,'PC_2'],
                    cell_type=seurat_control_KUL19$cell_type, state=seurat_control_KUL19$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/KUL19_samplings_reductions.csv')

# 7.7. Individual KUL21
seurat_control_KUL21 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$individual=='KUL21'])
invisible(gc())
# 7.7.1. Run PCA
seurat_control_KUL21 = Seurat::RunPCA(seurat_control_KUL21, features = Seurat::VariableFeatures(seurat_control_KUL21))
invisible(gc())
# 7.7.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_KUL21, reduction='pca', group.by='cell_type', cols=ct_colors) /
  Seurat::DimPlot(seurat_control_KUL21, reduction='pca', group.by='state', cols=state_colors)
# 7.7.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_KUL21@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_KUL21@reductions$pca@cell.embeddings[,'PC_2'],
                    cell_type=seurat_control_KUL21$cell_type, state=seurat_control_KUL21$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/KUL21_samplings_reductions.csv')

# 7.8. Individual SMC01
seurat_control_SMC01 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$individual=='SMC01'])
invisible(gc())
# 7.8.1. Run PCA
seurat_control_SMC01 = Seurat::RunPCA(seurat_control_SMC01, features = Seurat::VariableFeatures(seurat_control_SMC01))
invisible(gc())
# 7.8.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_SMC01, reduction='pca', group.by='cell_type', cols=ct_colors) /
  Seurat::DimPlot(seurat_control_SMC01, reduction='pca', group.by='state', cols=state_colors)
# 7.8.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_SMC01@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_SMC01@reductions$pca@cell.embeddings[,'PC_2'],
                    cell_type=seurat_control_SMC01$cell_type, state=seurat_control_SMC01$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/SMC01_samplings_reductions.csv')

# 7.9. Individual SMC04
seurat_control_SMC04 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$individual=='SMC04'])
invisible(gc())
# 7.9.1. Run PCA
seurat_control_SMC04 = Seurat::RunPCA(seurat_control_SMC04, features = Seurat::VariableFeatures(seurat_control_SMC04))
invisible(gc())
# 7.9.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_SMC04, reduction='pca', group.by='cell_type', cols=ct_colors) /
  Seurat::DimPlot(seurat_control_SMC04, reduction='pca', group.by='state', cols=state_colors)
# 7.9.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_SMC04@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_SMC04@reductions$pca@cell.embeddings[,'PC_2'],
                    cell_type=seurat_control_SMC04$cell_type, state=seurat_control_SMC04$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/SMC04_samplings_reductions.csv')

# 7.10. Individual SMC06
seurat_control_SMC06 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$individual=='SMC06'])
invisible(gc())
# 7.10.1. Run PCA
seurat_control_SMC06 = Seurat::RunPCA(seurat_control_SMC06, features = Seurat::VariableFeatures(seurat_control_SMC06))
invisible(gc())
# 7.10.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_SMC06, reduction='pca', group.by='cell_type', cols=ct_colors) /
  Seurat::DimPlot(seurat_control_SMC06, reduction='pca', group.by='state', cols=state_colors)
# 7.10.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_SMC06@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_SMC06@reductions$pca@cell.embeddings[,'PC_2'],
                    cell_type=seurat_control_SMC06$cell_type, state=seurat_control_SMC06$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/SMC06_samplings_reductions.csv')

# 7.11. Individual SMC07
seurat_control_SMC07 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$individual=='SMC07'])
invisible(gc())
# 7.11.1. Run PCA
seurat_control_SMC07 = Seurat::RunPCA(seurat_control_SMC07, features = Seurat::VariableFeatures(seurat_control_SMC07))
invisible(gc())
# 7.11.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_SMC07, reduction='pca', group.by='cell_type', cols=ct_colors) /
  Seurat::DimPlot(seurat_control_SMC07, reduction='pca', group.by='state', cols=state_colors)
# 7.11.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_SMC07@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_SMC07@reductions$pca@cell.embeddings[,'PC_2'],
                    cell_type=seurat_control_SMC07$cell_type, state=seurat_control_SMC07$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/SMC07_samplings_reductions.csv')

# 7.12. Individual SMC08
seurat_control_SMC08 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$individual=='SMC08'])
invisible(gc())
# 7.12.1. Run PCA
seurat_control_SMC08 = Seurat::RunPCA(seurat_control_SMC08, features = Seurat::VariableFeatures(seurat_control_SMC08))
invisible(gc())
# 7.12.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_SMC08, reduction='pca', group.by='cell_type', cols=ct_colors) /
  Seurat::DimPlot(seurat_control_SMC08, reduction='pca', group.by='state', cols=state_colors)
# 7.12.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_SMC08@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_SMC08@reductions$pca@cell.embeddings[,'PC_2'],
                    cell_type=seurat_control_SMC08$cell_type, state=seurat_control_SMC08$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/SMC08_samplings_reductions.csv')

# 7.13. Individual SMC10
seurat_control_SMC10 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$individual=='SMC10'])
invisible(gc())
# 7.13.1. Run PCA
seurat_control_SMC10 = Seurat::RunPCA(seurat_control_SMC10, features = Seurat::VariableFeatures(seurat_control_SMC10))
invisible(gc())
# 7.13.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_SMC10, reduction='pca', group.by='cell_type', cols=ct_colors) /
  Seurat::DimPlot(seurat_control_SMC10, reduction='pca', group.by='state', cols=state_colors)
# 7.13.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_SMC10@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_SMC10@reductions$pca@cell.embeddings[,'PC_2'],
                    cell_type=seurat_control_SMC10$cell_type, state=seurat_control_SMC10$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/SMC10_samplings_reductions.csv')


# 8. PCA for each cell-type separately
# 8.1. Cytotoxic CD8 Tcells
seurat_control_cyto = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$cell_type=='Cytotoxic CD8 Tcells'])
invisible(gc())
# 8.1.1. Run PCA
seurat_control_cyto = Seurat::RunPCA(seurat_control_cyto, features = Seurat::VariableFeatures(seurat_control_cyto))
invisible(gc())
# 8.1.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_cyto, reduction='pca', group.by='state', cols=state_colors)/
Seurat::DimPlot(seurat_control_cyto, reduction='pca', group.by='individual', cols=individual_colors)
# 8.1.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_cyto@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_cyto@reductions$pca@cell.embeddings[,'PC_2'],
                    individual=seurat_control_cyto$individual, state=seurat_control_cyto$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/cyto_samplings_reductions.csv')

# 8.2. Follicular CD4 Tcells
seurat_control_foll = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$cell_type=='Follicular CD4 Tcells'])
invisible(gc())
# 8.2.1. Run PCA
seurat_control_foll = Seurat::RunPCA(seurat_control_foll, features = Seurat::VariableFeatures(seurat_control_foll))
invisible(gc())
# 8.2.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_foll, reduction='pca', group.by='state', cols=state_colors)/
  Seurat::DimPlot(seurat_control_foll, reduction='pca', group.by='individual', cols=individual_colors)
# 8.2.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_foll@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_foll@reductions$pca@cell.embeddings[,'PC_2'],
                    individual=seurat_control_foll$individual, state=seurat_control_foll$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/foll_samplings_reductions.csv')

# 8.3. IL17+ CD4 Tcells
seurat_control_il17 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$cell_type=='IL17+ CD4 Tcells'])
invisible(gc())
# 8.3.1. Run PCA
seurat_control_il17 = Seurat::RunPCA(seurat_control_il17, features = Seurat::VariableFeatures(seurat_control_il17))
invisible(gc())
# 8.3.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_il17, reduction='pca', group.by='state', cols=state_colors)/
  Seurat::DimPlot(seurat_control_il17, reduction='pca', group.by='individual', cols=individual_colors)
# 8.3.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_il17@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_il17@reductions$pca@cell.embeddings[,'PC_2'],
                    individual=seurat_control_il17$individual, state=seurat_control_il17$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/il17_samplings_reductions.csv')

# 8.4. Memory CD4 Tcells
seurat_control_memcd4 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$cell_type=='Memory CD4 Tcells'])
invisible(gc())
# 8.4.1. Run PCA
seurat_control_memcd4 = Seurat::RunPCA(seurat_control_memcd4, features = Seurat::VariableFeatures(seurat_control_memcd4))
invisible(gc())
# 8.4.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_memcd4, reduction='pca', group.by='state', cols=state_colors)/
  Seurat::DimPlot(seurat_control_memcd4, reduction='pca', group.by='individual', cols=individual_colors)
# 8.4.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_memcd4@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_memcd4@reductions$pca@cell.embeddings[,'PC_2'],
                    individual=seurat_control_memcd4$individual, state=seurat_control_memcd4$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/memcd4_samplings_reductions.csv')

# 8.5. Memory CD8 Tcells
seurat_control_memcd8 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$cell_type=='Memory CD8 Tcells'])
invisible(gc())
# 8.5.1. Run PCA
seurat_control_memcd8 = Seurat::RunPCA(seurat_control_memcd8, features = Seurat::VariableFeatures(seurat_control_memcd8))
invisible(gc())
# 8.5.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_memcd8, reduction='pca', group.by='state', cols=state_colors)/
  Seurat::DimPlot(seurat_control_memcd8, reduction='pca', group.by='individual', cols=individual_colors)
# 8.5.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_memcd8@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_memcd8@reductions$pca@cell.embeddings[,'PC_2'],
                    individual=seurat_control_memcd8$individual, state=seurat_control_memcd8$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/memcd8_samplings_reductions.csv')

# 8.6. Naive CD4 Tcells
seurat_control_naivecd4 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$cell_type=='Naive CD4 Tcells'])
invisible(gc())
# 8.6.1. Run PCA
seurat_control_naivecd4 = Seurat::RunPCA(seurat_control_naivecd4, features = Seurat::VariableFeatures(seurat_control_naivecd4))
invisible(gc())
# 8.6.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_naivecd4, reduction='pca', group.by='state', cols=state_colors)/
  Seurat::DimPlot(seurat_control_naivecd4, reduction='pca', group.by='individual', cols=individual_colors)
# 8.6.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_naivecd4@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_naivecd4@reductions$pca@cell.embeddings[,'PC_2'],
                    individual=seurat_control_naivecd4$individual, state=seurat_control_naivecd4$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/naivecd4_samplings_reductions.csv')

# 8.7. Naive CD8 Tcells
seurat_control_naivecd8 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$cell_type=='Naive CD8 Tcells'])
invisible(gc())
# 8.7.1. Run PCA
seurat_control_naivecd8 = Seurat::RunPCA(seurat_control_naivecd8, features = Seurat::VariableFeatures(seurat_control_naivecd8))
invisible(gc())
# 8.7.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_naivecd8, reduction='pca', group.by='state', cols=state_colors)/
  Seurat::DimPlot(seurat_control_naivecd8, reduction='pca', group.by='individual', cols=individual_colors)
# 8.7.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_naivecd8@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_naivecd8@reductions$pca@cell.embeddings[,'PC_2'],
                    individual=seurat_control_naivecd8$individual, state=seurat_control_naivecd8$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/naivecd8_samplings_reductions.csv')

# 8.8. Proliferative CD4 Tcells
seurat_control_prolcd4 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$cell_type=='Proliferative CD4 Tcells'])
invisible(gc())
# 8.8.1. Run PCA
seurat_control_prolcd4 = Seurat::RunPCA(seurat_control_prolcd4, features = Seurat::VariableFeatures(seurat_control_prolcd4))
invisible(gc())
# 8.8.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_prolcd4, reduction='pca', group.by='state', cols=state_colors)/
  Seurat::DimPlot(seurat_control_prolcd4, reduction='pca', group.by='individual', cols=individual_colors)
# 8.8.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_prolcd4@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_prolcd4@reductions$pca@cell.embeddings[,'PC_2'],
                    individual=seurat_control_prolcd4$individual, state=seurat_control_prolcd4$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/prolcd4_samplings_reductions.csv')

# 8.9. Proliferative CD8 Tcells
seurat_control_prolcd8 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$cell_type=='Proliferative CD8 Tcells'])
invisible(gc())
# 8.9.1. Run PCA
seurat_control_prolcd8 = Seurat::RunPCA(seurat_control_prolcd8, features = Seurat::VariableFeatures(seurat_control_prolcd8))
invisible(gc())
# 8.9.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_prolcd8, reduction='pca', group.by='state', cols=state_colors)/
  Seurat::DimPlot(seurat_control_prolcd8, reduction='pca', group.by='individual', cols=individual_colors)
# 8.9.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_prolcd8@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_prolcd8@reductions$pca@cell.embeddings[,'PC_2'],
                    individual=seurat_control_prolcd8$individual, state=seurat_control_prolcd8$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/prolcd8_samplings_reductions.csv')

# 8.9. Regulatory CD4 Tcells
seurat_control_regcd4 = subset(seurat_control, cells=rownames(seurat_control[[]])[seurat_control$cell_type=='Regulatory CD4 Tcells'])
invisible(gc())
# 8.9.1. Run PCA
seurat_control_regcd4 = Seurat::RunPCA(seurat_control_regcd4, features = Seurat::VariableFeatures(seurat_control_regcd4))
invisible(gc())
# 8.9.2. PCA plots coloured by cell-types, state and individuals
Seurat::DimPlot(seurat_control_regcd4, reduction='pca', group.by='state', cols=state_colors)/
  Seurat::DimPlot(seurat_control_regcd4, reduction='pca', group.by='individual', cols=individual_colors)
# 8.9.3. Save dataframe with PCA information
infoDF = data.frame(PC1=seurat_control_regcd4@reductions$pca@cell.embeddings[,'PC_1'],
                    PC2=seurat_control_regcd4@reductions$pca@cell.embeddings[,'PC_2'],
                    individual=seurat_control_regcd4$individual, state=seurat_control_regcd4$state)
write.csv(infoDF, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/regcd4_samplings_reductions.csv')





# -----
# --- Density plot of relevant paths for nive and proliferative cells
# -----

important_rxns = jsonlite::read_json('./GENERAL/utility_data/important_reactions_Tcells.json',
                                     simplifyVector=TRUE)
naive_prol = seurat_control$cell_type%in%c('Naive CD4 Tcells', 'Naive CD8 Tcells', 'Proliferative CD4 Tcells', 'Proliferative CD8 Tcells')

# 1. Glycolysis:
glycolysis_rxns = c('MAR04394', 'MAR07747', 'MAR04381', 'MAR04301', 'MAR04379',
                    'MAR04375', 'MAR04391', 'MAR04373', 'MAR04371', 'MAR04372',
                    'MAR04368', 'MAR04370', 'MAR04365', 'MAR04363', 'MAR04378')
plts=list()
for(rxn in glycolysis_rxns){
  message(rxn)
  df = data.frame(reaction=seurat_control@assays$RNA@counts[rxn, naive_prol], cell_type=seurat_control$cell_type[naive_prol])
  plts[[rxn]]=ggplot2::ggplot(df, ggplot2::aes(reaction, y=..scaled.., colour=cell_type)) +
    ggplot2::geom_density(adjust=.1) + ggplot2::scale_colour_manual(values=ct_colors) +
    ggplot2::theme(legend.position = 'none') + ggplot2::xlab(rxn) + ggplot2::ylab('')
}
patchwork::wrap_plots(plts, ncol=3)

# 2. Pyruvate to lactate vs pyruvte to acetyl-CoA:
pyrConsumption_rxns = c('MAR04388', 'MAR04137')
plts=list()
for(rxn in pyrConsumption_rxns){
  message(rxn)
  df = data.frame(reaction=seurat_control@assays$RNA@counts[rxn, naive_prol], cell_type=seurat_control$cell_type[naive_prol])
  plts[[rxn]]=ggplot2::ggplot(df, ggplot2::aes(reaction, y=..scaled.., colour=cell_type)) +
    ggplot2::geom_density(adjust=.1) + ggplot2::scale_colour_manual(values=ct_colors) +
    ggplot2::theme(legend.position = 'none') + ggplot2::xlab(rxn) + ggplot2::ylab('')
}
patchwork::wrap_plots(plts, ncol=1)

# 3. OXPHOS:
oxphos_rxns = c('MAR06921', 'MAR06911', 'MAR06918', 'MAR06914', 'MAR13081', 'MAR06916')
plts=list()
for(rxn in oxphos_rxns){
  message(rxn)
  df = data.frame(reaction=seurat_control@assays$RNA@counts[rxn, naive_prol], cell_type=seurat_control$cell_type[naive_prol])
  plts[[rxn]]=ggplot2::ggplot(df, ggplot2::aes(reaction, y=..scaled.., colour=cell_type)) +
    ggplot2::geom_density(adjust=.1) + ggplot2::scale_colour_manual(values=ct_colors) +
    ggplot2::theme(legend.position = 'none') + ggplot2::xlab(rxn) + ggplot2::ylab('')
}
patchwork::wrap_plots(plts, ncol=3)

# 4. TCA:
tca_rxns = c('MAR04145', 'MAR04458', 'MAR04589', 'MAR03957', 'MAR03958', 'MAR04588', 'MAR04112',
             'MAR05297', 'MAR04209', 'MAR06413', 'MAR06411', 'MAR06414', 'MAR03787', 'MAR04147',
             'MAR04152', 'MAR04652', 'MAR08743', 'MAR04410', 'MAR04141', 'MAR04145')
plts=list()
for(rxn in tca_rxns){
  message(rxn)
  df = data.frame(reaction=seurat_control@assays$RNA@counts[rxn, naive_prol], cell_type=seurat_control$cell_type[naive_prol])
  plts[[rxn]]=ggplot2::ggplot(df, ggplot2::aes(reaction, y=..scaled.., colour=cell_type)) +
    ggplot2::geom_density(adjust=.1) + ggplot2::scale_colour_manual(values=ct_colors) +
    ggplot2::theme(legend.position = 'none') + ggplot2::xlab(rxn) + ggplot2::ylab('')
}
patchwork::wrap_plots(plts, ncol=3)

# 5. FAS
#col_annot = ComplexHeatmap::columnAnnotation(Cell_Type=seurat_control$cell_type[cells_use],
#                                             State=seurat_control$state[cells_use],
#                                             Individual=seurat_control$individual[cells_use],
#                                             col=list(Cell_Type=ct_colors, State=state_colors,
#                                                      Individual=individual_colors))
#FAS_reactions = important_rxns$subsystems$`Fatty acid biosynthesis`[as.vector(important_rxns$subsystems$`Fatty acid biosynthesis`)%in%rownames(seurat_control@assays$RNA@counts)]
#ComplexHeatmap::Heatmap(as.matrix(seurat_control@assays$RNA@counts[FAS_reactions, cells_use]),
#                       bottom_annotation=col_annot,
#                        col=circlize::colorRamp2(c(-1,0,1), c('blue', 'white', 'red')),
#                        column_split=seurat_control$cell_type[cells_use], column_title=NULL,
#                        cluster_columns=FALSE, cluster_rows=TRUE,
#                        show_column_names=FALSE, show_row_names=FALSE)
#

# 6. Glutaminolysis
glutaminolysis_reactions = c('MAR03892', 'MAR03827', 'MAR04109', 'MAR03802', 'MAR03804')
plts=list()
for(rxn in glutaminolysis_reactions){
  message(rxn)
  df = data.frame(reaction=seurat_control@assays$RNA@counts[rxn, naive_prol], cell_type=seurat_control$cell_type[naive_prol])
  plts[[rxn]]=ggplot2::ggplot(df, ggplot2::aes(reaction, y=..scaled.., colour=cell_type)) +
    ggplot2::geom_density(adjust=.1) + ggplot2::scale_colour_manual(values=ct_colors) +
    ggplot2::theme(legend.position = 'none') + ggplot2::xlab(rxn) + ggplot2::ylab('')
}
patchwork::wrap_plots(plts, ncol=3)

# 7. BCAA catabolism
plts=list()
for(rxn in important_rxns$subsystems$`Branched-chain amino acids (BCAA) catabolism`){
  message(rxn)
  df = data.frame(reaction=seurat_control@assays$RNA@counts[rxn, naive_prol], cell_type=seurat_control$cell_type[naive_prol])
  plts[[rxn]]=ggplot2::ggplot(df, ggplot2::aes(reaction, y=..scaled.., colour=cell_type)) +
    ggplot2::geom_density(adjust=.1) + ggplot2::scale_colour_manual(values=ct_colors) +
    ggplot2::theme(legend.position = 'none') + ggplot2::xlab(rxn) + ggplot2::ylab('')
}
patchwork::wrap_plots(plts, ncol=3)

# 8. Glucose in PPP
plts=list()
glucose_in_PPP = c('MAR04306', 'MAR08971')
for(rxn in glucose_in_PPP){
  message(rxn)
  df = data.frame(reaction=seurat_control@assays$RNA@counts[rxn, ], cell_type=seurat_control$cell_type)
  plts[[rxn]]=ggplot2::ggplot(df, ggplot2::aes(reaction, y=..scaled.., colour=cell_type)) +
    ggplot2::geom_density(adjust=.1) + ggplot2::scale_colour_manual(values=ct_colors) +
    ggplot2::theme(legend.position = 'none') + ggplot2::xlab(rxn) + ggplot2::ylab('')
}
patchwork::wrap_plots(plts, ncol=2)

# 9. Creation of precursors for purines, pyrimidines
plts=list()
purine_pyrimidine_precursors = c('MAR04038', 'MAR04577')
for(rxn in purine_pyrimidine_precursors){
  message(rxn)
  df = data.frame(reaction=seurat_control@assays$RNA@counts[rxn, ], cell_type=seurat_control$cell_type)
  plts[[rxn]]=ggplot2::ggplot(df, ggplot2::aes(reaction, y=..scaled.., colour=cell_type)) +
    ggplot2::geom_density(adjust=.1) + ggplot2::scale_colour_manual(values=ct_colors) +
    ggplot2::theme(legend.position = 'none') + ggplot2::xlab(rxn) + ggplot2::ylab('')
}
patchwork::wrap_plots(plts, ncol=2)

# 10. Protein synthesis
df = data.frame(reaction=seurat_control@assays$RNA@counts['MAR13078', ], cell_type=seurat_control$cell_type)
ggplot2::ggplot(df, ggplot2::aes(reaction, y=..scaled.., colour=cell_type)) +
  ggplot2::geom_density(adjust=.1) + ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::theme(legend.position = 'none') + ggplot2::xlab('MAR13078') + ggplot2::ylab('')

# 11. GLUT1
plts=list()
for(rxn in important_rxns$uptakes$GLUT1){
  message(rxn)
  df = data.frame(reaction=seurat_control@assays$RNA@counts[rxn, ], cell_type=seurat_control$cell_type)
  plts[[rxn]]=ggplot2::ggplot(df, ggplot2::aes(reaction, y=..scaled.., colour=cell_type)) +
    ggplot2::geom_density(adjust=.1) + ggplot2::scale_colour_manual(values=ct_colors) +
    ggplot2::theme(legend.position = 'none') + ggplot2::xlab(rxn) + ggplot2::ylab('')
}
patchwork::wrap_plots(plts, ncol=3)

# 12. HK2
plts=list()
for(rxn in important_rxns$HK2){
  message(rxn)
  df = data.frame(reaction=seurat_control@assays$RNA@counts[rxn, ], cell_type=seurat_control$cell_type)
  plts[[rxn]]=ggplot2::ggplot(df, ggplot2::aes(reaction, y=..scaled.., colour=cell_type)) +
    ggplot2::geom_density(adjust=.1) + ggplot2::scale_colour_manual(values=ct_colors) +
    ggplot2::theme(legend.position = 'none') + ggplot2::xlab(rxn) + ggplot2::ylab('')
}
patchwork::wrap_plots(plts, ncol=3)

# 13. FPK1 and FPK2
df = data.frame(reaction=seurat_control@assays$RNA@counts['MAR04379', ], cell_type=seurat_control$cell_type)
ggplot2::ggplot(df, ggplot2::aes(reaction, y=..scaled.., colour=cell_type)) +
  ggplot2::geom_density(adjust=.1) + ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::theme(legend.position = 'none') + ggplot2::xlab('MAR04379') + ggplot2::ylab('')

# 14. MCT1
df = data.frame(reaction=seurat_control@assays$RNA@counts['MAR01517', ], cell_type=seurat_control$cell_type)
ggplot2::ggplot(df, ggplot2::aes(reaction, y=..scaled.., colour=cell_type)) +
  ggplot2::geom_density(adjust=.1) + ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::theme(legend.position = 'none') + ggplot2::xlab('MAR01517') + ggplot2::ylab('')

# 15. LDHA
plts=list()
for(rxn in important_rxns$LDHA){
  message(rxn)
  df = data.frame(reaction=seurat_control@assays$RNA@counts[rxn, ], cell_type=seurat_control$cell_type)
  plts[[rxn]]=ggplot2::ggplot(df, ggplot2::aes(reaction, y=..scaled.., colour=cell_type)) +
    ggplot2::geom_density(adjust=.1) + ggplot2::scale_colour_manual(values=ct_colors) +
    ggplot2::theme(legend.position = 'none') + ggplot2::xlab(rxn) + ggplot2::ylab('')
}
patchwork::wrap_plots(plts, ncol=3)





pcor = cor(as.matrix(Matrix::t(seurat_control@assays$RNA@counts)), method='pearson')
write.csv(pcor, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/correlation.csv')

x=which(pcor>0.9, arr.ind=T)
x[x[,1]<x[,2],]

x2=which(pcor<(-0.9), arr.ind=T)
x2[x2[,1]<x2[,2],]

ComplexHeatmap::Heatmap(pcor, cluster_columns=FALSE, cluster_rows=FALSE,
                        show_column_names=FALSE, show_row_names=FALSE)






# -----
# --- Compare important reactions between cell-types + state (e.g. naiveCD4_normal)
# -----

sampling_control_data = Matrix::readMM('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/sampling_control_data.mtx')
CRCatlas_samplingNMControl_meta = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/sampling_control_metadata.csv',
                                           row.names=1, stringsAsFactors=TRUE)

important_rxns = jsonlite::read_json('./GENERAL/utility_data/important_reactions_Tcells.json',
                                     simplifyVector=TRUE)

# Biomass:
regs = CRCatlas_samplingNMControl_meta$cell_type=='Regulatory CD4 Tcells'
density_rxn_StateCT(sampling_control_data[regs,], CRCatlas_samplingNMControl_meta[regs,],
                    important_rxns$biomass, colour_by=NULL, cts_order=NULL,
                    state_order=NULL)

# SLC7A5:
density_sumrxns_StateCT(sampling_control_data, CRCatlas_samplingNMControl_meta,
                        important_rxns$uptakes$SLC7A5, 'SLC7A5',
                        colour_by=NULL, cts_order=NULL, state_order=NULL,
                        abs=T)

# LDHA:
density_sumrxns_StateCT(sampling_control_data, CRCatlas_samplingNMControl_meta,
                        important_rxns$LDHA, 'LDHA', colour_by=NULL, cts_order=NULL,
                        state_order=NULL)

# Glycolysis pathway scores:
glycolysis_scores = glycolysis_score(sampling_control_data)
density_pathscores_StateCT(glycolysis_scores, CRCatlas_samplingNMControl_meta,
                           'Glycolysis', colour_by=NULL, cts_order=NULL,
                           state_order=NULL)

# Lactate/TCA ratio from pyruvate:
lac_tca = pyruvate_to_lactate_vs_TCA(CRCatlas_samplingNMControl_matrix)
lac_tca2 = lac_tca
lac_tca2[is.na(lac_tca2)] = 0
boxplot_pathscores_StateCT(lac_tca, CRCatlas_samplingNMControl_meta,
                           'Ratio lactate vs TCA', 'state',
                           colour_by='cell_type')






# -----
# --- Pre-Process data for differential and PCA analyses
# -----

# Read matrix:
sampling_control_data = Matrix::readMM('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/sampling_control_data.mtx')

# Transpose matrix
sampling_control_data = Matrix::t(sampling_control_data)
# Now it is reactions vs samplings, i.e., features vs samples.

# Remove reactions with non-zero flux in less than 100 samplings
nonzero = sampling_control_data > 0
keep_rxns = Matrix::rowSums(nonzero) >= 100
sampling_control_data = sampling_control_data[keep_rxns, ]
remove(nonzero)
invisible(gc())

# Normalize samplings to sum to 1 and find highly variable reactions (followed
# Seurats' code)
# Normalization by sum
norm_factors = Matrix::colSums(sampling_control_data)
sampling_control_dataNorm = Matrix::Matrix(t(apply(sampling_control_data, 1,
                                                   function(x) x/norm_factors)),
                                           sparse=T)
invisible(gc())
# Calculate reactions means, variances
rxn_means = Matrix::rowMeans(sampling_control_dataNorm)
rxn_vars = Seurat:::SparseRowVar2(mat=sampling_control_dataNorm,
                                  mu=rxn_means, display_progress=T)
# Standardized variance values must be maximum clip.max.
# Those bigger than clip.max will be set to clip.max
clip.max = sqrt(ncol(sampling_control_dataNorm))
# Reactions that are not constant:
rxns_notConst = rxn_vars > 0
# Fit a regression line (loess)
data.loess.info = data.frame(mean=rxn_means, variance=rxn_vars,
                        variance.expected=0, variance.standardized=0)
loess.fit = loess(formula=log10(rxn_vars) ~ log10(rxn_means),
                  data=data.loess.info, span=0.3)
# Save fit results in data.loess.info
data.loess.info$variance.expected[rxns_notConst] <- 10 ^ loess.fit$fitted
data.loess.info$variance.standardized = 
  Seurat:::SparseRowVarStd(mat=sampling_control_dataNorm, mu=data.loess.info$mean,
                           sd=sqrt(data.loess.info$variance.expected),
                           vmax=clip.max, display_progress=T)
# Get top 2 000 most variable reactions
data.loess.info = data.loess.info[which(data.loess.info[, 1, drop=TRUE]!=0), ]
data.loess.info = data.loess.info[order(data.loess.info$variance.standardized,
                                        decreasing=TRUE), , drop=FALSE]
top_rxns = head(rownames(data.loess.info), n=2000)

# Filter raw matrix to have only the highly variable reactions
sampling_control_data = sampling_control_data[top_rxns, ]
invisible(gc())

# Scale raw matrix to for mean = 0 and variance = 1
sampling_control_dataScaled = scale(sampling_control_data)

#sampling_control_AllDataScaled = scale(Matrix::t(sampling_control_data_original))

# BATCH CORRECTION FOR DATASETS? WHERE IN THE PROCESSING PIPELINE?

# Save scaled data
Matrix::writeMM(sampling_control_dataScaled,
                './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/sampling_control_dataScaled.mtx')





# -----
# --- UMAP
# -----

# PCA on scaled data, with weighting of feature embeddings by the variance of each PC
pca.results = irlba::irlba(A=t(sampling_control_dataScaled), nv=50)
rxn_embeddings = pca.results$u %*% diag(pca.results$d)

pca_plot_df = cbind(rxn_embeddings[,1:2], CRCatlas_samplingNMControl_meta)
colnames(pca_plot_df)[1:2] = c('PCA1', 'PCA2')
plt = ggplot2::ggplot(pca_plot_df, ggplot2::aes(PCA1, PCA2))
plt + ggplot2::geom_point(ggplot2::aes_string(colour='cell_type'))
plt + ggplot2::geom_point(ggplot2::aes_string(colour='medium'))
plt + ggplot2::geom_point(ggplot2::aes_string(colour='sample'))

# Where components start to elbow
#stdev = pca.results$d/sqrt(max(1, ncol(sampling_control_dataScaled) - 1))
#pct = stdev / sum(stdev) * 100
#cumu = cumsum(pct)
#co1 = which(cumu > 90 & pct < 5)[1]
#co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
#elbow = min(co1, co2)
#plot_df = data.frame(pct=pct, cumu=cumu, rank=1:length(pct))
#ggplot2::ggplot(plot_df, ggplot2::aes(cumu, pct, label = rank, color = rank > elbow)) + 
#  ggplot2::geom_text() + 
#  ggplot2::geom_vline(xintercept = 90, color = "grey") + 
#  ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
#  ggplot2::theme_bw()
#
# UMAP
#CRCatlas_samplingNMControl_umap = uwot::umap(rxn_embeddings[,1:elbow],
#                                             n_neighbors = 30, init='spectral',
#                                             metric='cosine', verbose = TRUE)
#
# Plots of result
#colnames(CRCatlas_samplingNMControl_umap) = c('UMAP_1', 'UMAP_2')
#df_umapPlot = cbind(CRCatlas_samplingNMControl_umap,
#                    CRCatlas_samplingNMControl_meta)
#plot_umap(df_umapPlot, 'cell_type', 'medium')
#plot_umap(df_umapPlot, 'medium')
#plot_umap(df_umapPlot, 'sample')




# -----
# --- Differential Flux Analysis
# -----

# Reactions whose fluxes are significantly different between a cell-type and 
# all others, for each medium

ct_results = list()
for(medium in c('blood', 'plasmax')){
  message(medium)
  ct_results[[medium]] = c()
  for(ct in unique(CRCatlas_samplingNMControl_meta$cell_type)){
    message('-',ct)
    medium_samplings = CRCatlas_samplingNMControl_meta$medium==medium
    ct_medium_samplings = CRCatlas_samplingNMControl_meta$cell_type==ct & medium_samplings
    other_medium_samplings = CRCatlas_samplingNMControl_meta$cell_type!=ct & medium_samplings
    # - Only consider reactions whose |log2(fold change)| is higher than 1
    ct_rxns_means = rowMeans(sampling_control_dataScaled[,ct_medium_samplings])
    other_rxns_means = rowMeans(sampling_control_dataScaled[,other_medium_samplings])
    fc = ct_rxns_means - other_rxns_means
    names(fc) = row.names(sampling_control_dataScaled)
    rxns_toTest = names(fc)[fc>1]
    # - For those reactions, find the ones that are significantly different
    test_meta1 = CRCatlas_samplingNMControl_meta$cell_type
    test_meta1 = data.frame(cell_type=test_meta1[medium_samplings],
                            row.names=colnames(sampling_control_dataScaled)[medium_samplings])
    test_meta1$cell_type[test_meta1$cell_type!=ct] = 'other'
    p_val = sapply(rxns_toTest,
      FUN = function(x) {
        return(wilcox.test(sampling_control_dataScaled[x, medium_samplings] ~ test_meta1[, "cell_type"])$p.value)
      }
    )
    pval_adjust = p.adjust(p_val, method='fdr')
    res_ct = data.frame(reaction=rxns_toTest, pval_adjust,
                        fold_change=fc[rxns_toTest], pval=p_val,
                        cell_type=ct, row.names=paste(rxns_toTest, ct, sep='_'))
    ct_results[[medium]] = rbind(ct_results[[medium]], res_ct)
  }
}

# ----------------------------------------------------------------------------
# Visualize differential reactions
rxn = 'MAR07881'
ct = 'foll'
ct_samplings = CRCatlas_samplingNMControl_meta$cell_type==ct
df = data.frame(reaction=sampling_control_dataScaled[rxn, ],
                cell_type=CRCatlas_samplingNMControl_meta$cell_type)
df$cell_type = factor(df$cell_type, levels=unique(df$cell_type))
ggplot2::ggplot(df, ggplot2::aes(reaction, colour=cell_type)) +
  ggplot2::geom_density() + 
  ggplot2::geom_vline(xintercept=mean(sampling_control_dataScaled[rxn,ct_samplings]), colour='red',
                      linetype='dashed') +
  ggplot2::geom_vline(xintercept=mean(sampling_control_dataScaled[rxn,!ct_samplings]), colour='gray',
                      linetype='dashed') +
  ggplot2::xlab(rxn)

# Ordered list of subsystems with most amount of significant reactions
ct='prolCD8'
significant_rxns = ct_results$reaction[ct_results$pval_adjust < 0.01 & ct_results$cell_type==ct]
sort(table(unlist(HumanGEM_rxns_subsystems[significant_rxns])))

all_significant_rxns = unique(ct_results$reaction)
subsystems_annotation = unlist(HumanGEM_rxns_subsystems[all_significant_rxns])
heatmapData = c()
for(ct in unique(CRCatlas_samplingNMControl_meta$cell_type)){
  ct_rxns_means = rowMeans(sampling_control_dataScaled[all_significant_rxns, CRCatlas_samplingNMControl_meta$cell_type==ct])
  other_rxns_means = rowMeans(sampling_control_dataScaled[all_significant_rxns, CRCatlas_samplingNMControl_meta$cell_type!=ct])
  fc = ct_rxns_means - other_rxns_means
  heatmapData = cbind(heatmapData, fc)
}
colnames(heatmapData) = unique(CRCatlas_samplingNMControl_meta$cell_type)
heatmapmatrix = as.matrix(heatmapData)
simple_subsystems_annotation = sapply(subsystems_annotation, function(x) strsplit(x, ' [(]')[[1]][1])

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
top10_represented_subsystems = names(sort(table(simple_subsystems_annotation), decreasing=T)[2:10])
row_annot=ComplexHeatmap::rowAnnotation(subsystem=simple_subsystems_annotation[simple_subsystems_annotation%in%top10_represented_subsystems])
ComplexHeatmap::Heatmap(heatmapmatrix[simple_subsystems_annotation%in%top10_represented_subsystems,],
                        right_annotation = row_annot, row_split = simple_subsystems_annotation[simple_subsystems_annotation%in%top10_represented_subsystems],
                        cluster_columns=FALSE, cluster_rows=TRUE,
                        show_column_names=TRUE, show_row_names=FALSE)


#rand_samps = sort(sample(1:10000, 5000))
#col_annot=ComplexHeatmap::HeatmapAnnotation(cell_type=CRCatlas_samplingNMControl_meta$cell_type[rand_samps])
#ht1 = ComplexHeatmap::Heatmap(as.matrix(sampling_control_data[res$pval_adjust<0.01, rand_samps]),
#                              cluster_columns=FALSE, cluster_rows=TRUE,
#                              column_split=CRCatlas_samplingNMControl_meta$cell_type[rand_samps],
#                              show_column_names=FALSE, show_row_names=FALSE,
#                              top_annotation=col_annot)
#ht1 = ComplexHeatmap::draw(ht1)
#rowOrder = ComplexHeatmap::row_order(ht1)
#
#col_annot=ComplexHeatmap::HeatmapAnnotation(cell_type=CRCatlas_samplingNMControl_meta$cell_type[rand_samps])
#ComplexHeatmap::Heatmap(as.matrix(sampling_control_dataScaled[res$pval_adjust<0.01, rand_samps]),
#                        cluster_columns=FALSE, cluster_rows=TRUE, row_order=rowOrder,
#                        column_split=CRCatlas_samplingNMControl_meta$cell_type[rand_samps],
#                        show_column_names=FALSE, show_row_names=FALSE,
#                        top_annotation=col_annot)
#
#
#
#smaller_than250 = sampling_control_data < 100 & sampling_control_data > -100
#rxns_smaller_than250 = Matrix::rowSums(smaller_than250) >= 9000
#col_annot=ComplexHeatmap::HeatmapAnnotation(cell_type=CRCatlas_samplingNMControl_meta$cell_type[rand_samps])
#ComplexHeatmap::Heatmap(as.matrix(sampling_control_data[res$pval_adjust<0.01 & rxns_smaller_than250, rand_samps]),
#                        cluster_columns=FALSE, cluster_rows=TRUE,
#                        column_split=CRCatlas_samplingNMControl_meta$cell_type[rand_samps],
#                        show_column_names=FALSE, show_row_names=FALSE,
#                        top_annotation=col_annot)
# ----------------------------------------------------------------------------


# Reactions whose fluxes are significantly different between the two media,
# for each cell-type


