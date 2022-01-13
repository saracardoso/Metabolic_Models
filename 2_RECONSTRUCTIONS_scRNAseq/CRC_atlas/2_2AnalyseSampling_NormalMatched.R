code_dir = './GENERAL/code/R'
for(file in list.files(code_dir, full.names=TRUE)) source(file)



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
                                                              'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['31']]$control$scrEXT002 = c('Cytotoxic CD8 Tcells',
                                                              'Memory CD4 Tcells',
                                                              'Proliferative CD8 Tcells',
                                                              'Follicular CD4 Tcells')
CRCatlas_sampling$NormalMatched[['31']]$control$scrEXT003 = c('Naive CD4 Tcells')

# Patient 32
CRCatlas_sampling$NormalMatched[['32']] = list()
CRCatlas_sampling$NormalMatched[['32']]$control = list()
CRCatlas_sampling$NormalMatched[['32']]$control$scrEXT009 = c('Cytotoxic CD8 Tcells',
                                                              'IL17+ CD4 Tcells', 'Memory CD8 Tcells',
                                                              'Naive CD8 Tcells',
                                                              'Proliferative CD4 Tcells',
                                                              'Proliferative CD8 Tcells')
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
CRCatlas_sampling$NormalMatched[['33']]$control$scrEXT012 = c('Memory CD8 Tcells',
                                                              'Proliferative CD4 Tcells',
                                                              'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['33']]$control$scrEXT013 = c('Memory CD4 Tcells',
                                                              'Naive CD4 Tcells',
                                                              'Proliferative CD4 Tcells')
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
                                                                    'Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
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
CRCatlas_sampling$NormalMatched[['KUL21']]$control[['KUL21-B']] = c('Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells')
CRCatlas_sampling$NormalMatched[['KUL21']]$control[['KUL21-N']] = c('IL17+ CD4 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['KUL21']]$control[['KUL21-T']] = c('Memory CD8 Tcells',
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
                                                                    'Naive CD8 Tcells')
CRCatlas_sampling$NormalMatched[['SMC06']]$control[['SMC06-T']] = c('Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')

# Patient SMC07
CRCatlas_sampling$NormalMatched[['SMC07']] = list()
CRCatlas_sampling$NormalMatched[['SMC07']]$control = list()
CRCatlas_sampling$NormalMatched[['SMC07']]$control[['SMC07-N']] = c('Memory CD4 Tcells',
                                                                    'Naive CD4 Tcells')
CRCatlas_sampling$NormalMatched[['SMC07']]$control[['SMC07-T']] = c('Cytotoxic CD8 Tcells',
                                                                    'Follicular CD4 Tcells',
                                                                    'IL17+ CD4 Tcells',
                                                                    'Memory CD4 Tcells',
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
                                                                    'Proliferative CD8 Tcells')

# Save the list into a json file:
jsonlite::write_json(CRCatlas_sampling, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/CRCatlas_sampling.json')





# -----
# --- Read sampling data and metadata - NormalMatched + control
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
  for(medium_type in c('Blood_SMDB', 'Plasmax_serum')){
    # Add to list:
    mtx = Matrix::Matrix(as.matrix(read.csv(paste('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched',
                                                  individual, '1_sampling/control', sample, gsub(' ', '_', cell_type),
                                                  paste(medium_type, '.csv', sep=''), sep='/'), row.names=1)), sparse=T)
    # Change rownames:
    rownms = paste(sample, gsub(' ', '_', cell_type), medium_type, rownames(mtx), sep='_')
    rownames(mtx) = rownms
    
    # Add to matrix:
    CRCatlas_samplingNMControl_matrix = rbind(CRCatlas_samplingNMControl_matrix, mtx)
    
    # Save data until now:
    Matrix::writeMM(CRCatlas_samplingNMControl_matrix,
                    './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/sampling_control_data.mtx')
    
    # Samples metadata: [sample, individual, state, cell_type]
    n_samplings = dim(CRCatlas_samplingNMControl[[individual]][[cell_type]][[medium_type]])[1]
    samples = c(samples, rep(sample, n_samplings))
    individuals = c(individuals, rep(individual, n_samplings))
    states = c(states, rep(metadata[sample, 'Sample.Source'], n_samplings))
    cell_types = c(cell_types, rep(cell_type, n_samplings))
    media = c(media, rep(medium_type, n_samplings))
    row_names = c(row_names, rownms)
  }
}
CRCatlas_samplingNMControl_meta = data.frame(sample=samples, individual=individuals, state=states,
                                             cell_type=cell_types, medium=media, row.names=row_names)

# Save metadata:
write.csv(CRCatlas_samplingNMControl_meta,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/sampling_control_metadata.csv')



# -----
# --- UMAP
# -----

# Scale

# PCA

# Run UMAP on PCA
CRCatlas_samplingNMControl_umap = umap_sampling(CRCatlas_samplingNMControl) # this function will need to change!!

# Plots result
df_umapPlot = cbind(CRCatlas_samplingNMControl_umap, CRCatlas_samplingNMControl_meta)
plot_umap(df_umapPlot, '')



