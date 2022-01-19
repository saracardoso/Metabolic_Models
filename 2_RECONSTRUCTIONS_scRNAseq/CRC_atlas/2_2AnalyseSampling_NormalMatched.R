code_dir = './GENERAL/code/R'
for(file in list.files(code_dir, full.names=TRUE)) source(file)

HumanGEM_rxns_subsystems = jsonlite::read_json('./GENERAL/utility_data/reactions_subsystems_mapping.json',
                                               simplifyVector = TRUE)
HumanGEM_subsystems_rxns = jsonlite::read_json('./GENERAL/utility_data/subsystems_reactions_mapping.json',
                                               simplifyVector = TRUE)



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
  for(medium_type in c('Blood_SMDB', 'Plasmax_serum')){
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
                    './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/sampling_control_data.mtx')
    
    # Samples metadata: [sample, individual, state, cell_type]
    message('- Collecting metadata...')
    n_samplings = dim(mtx)[1]
    samples = c(samples, rep(sample, n_samplings))
    individuals = c(individuals, rep(individual, n_samplings))
    states = c(states, rep(metadata[sample, 'Sample.Source'], n_samplings))
    cell_types = c(cell_types, rep(cell_type, n_samplings))
    media = c(media, rep(medium_type, n_samplings))
    row_names = c(row_names, rownms)
    
    # Deleting temporary mtx:
    message('- Deleting temporary mtx...')
    remove(mtx)
    invisible(gc())
  }
  message('')
}
CRCatlas_samplingNMControl_meta = data.frame(sample=samples, individual=individuals, state=states,
                                             cell_type=cell_types, medium=media, row.names=row_names)

# Save metadata:
write.csv(CRCatlas_samplingNMControl_meta,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/sampling_control_metadata.csv')





# -----
# --- Pre-Process data for analysis
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


