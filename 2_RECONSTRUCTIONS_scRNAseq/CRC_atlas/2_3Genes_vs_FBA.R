ct_colors = c('#b30000', '#999999', '#cc79a7', '#b3c9e6', '#4e84c4', '#c3d7a4', '#52854c', '#fdb981',
              '#d16103', '#293352')
names(ct_colors) = c('Cytotoxic CD8 Tcells', 'Follicular CD4 Tcells', 'IL17+ CD4 Tcells',
                     'Memory CD4 Tcells', 'Memory CD8 Tcells', 'Naive CD4 Tcells', 'Naive CD8 Tcells',
                     'Proliferative CD4 Tcells', 'Proliferative CD8 Tcells',  'Regulatory CD4 Tcells')
state_colors = c('#52854c', '#b30000', '#ff6666', '#cc79a7')
names(state_colors) = c('Normal Matched', 'Tumour', 'Tumour Core', 'Tumour Border')

# ----------------------------------------------
# --- scRNAseq data used for models analysed ---
# ----------------------------------------------

# Read scRNAseq CRC atlas:
tcells = SeuratDisk::LoadH5Seurat('/home/scardoso/Documents/PhD/CRC_ATLAS/2_annotation/results_Tcells/datasets/Tcells_finalAnnots.h5Seurat')

# Get only individuals that were used in the models:
indivs_keep = rownames(tcells[[]])[!tcells$patient%in%c('36', '37', '38', 'KUL28', 'KUL30', 'KUL31',
                                                           'SMC02', 'SMC03', 'SMC05', 'SMC09', 'N13')]
tcells = subset(tcells, cells=indivs_keep)
invisible(gc())

# Get cell-type annotation equal to the models:
tcells[['models_annotation']] = tcells$Annotation_Level_4

# Get only the cells from the models that ended up being analysed:
models_created = jsonlite::read_json('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/CRCatlas_sampling.json',
                                     simplifyVector=TRUE)
cells_keep = c()
for(state in models_created){
  for(indiv in state){
    for(samp in names(indiv$control)){
      models = indiv$control[[samp]]
      cells_keep = c(cells_keep,
                     rownames(tcells[[]])[tcells$sample==samp & tcells$models_annotation%in%models])
    }
  }
}
tcells = subset(tcells, cells=cells_keep)
invisible(gc())



# --------------------------------------------------------------------
# --- scRNAseq data used for models analysed, only metabolic genes ---
# --------------------------------------------------------------------

humangem_genemapping = jsonlite::read_json('./GENERAL/utility_data/genes_mapping.json', simplifyVector=TRUE)
metabolic_genes = names(humangem_genemapping)
metabolic_genes_scrnaseq = metabolic_genes[metabolic_genes%in%rownames(tcells@assays$RNA@counts)]
tcells_metabolic = subset(tcells, features = metabolic_genes_scrnaseq)
invisible(gc())



# -----------------------------------------------------------------------
# --- pseudo-bulk data used for models analysed, only metabolic genes ---
# -----------------------------------------------------------------------

pseudo_bulk_tcells_allGenes = NULL
tcelltypes_metabolic = readRDS('./0Data/scRNAseq/CRC_atlas/expression_data.Rdata')
for(indiv in names(tcelltypes_metabolic)){
  if(indiv%in%names(models_created$NormalMatched)){
    for(samp in names(tcelltypes_metabolic[[indiv]])){
      message(samp)
      expression = tcelltypes_metabolic[[indiv]][[samp]]
      expression = expression[, colnames(expression)%in%models_created$NormalMatched[[indiv]]$control[[samp]]]
      expression = cbind(rownames(expression), expression)
      colnames(expression) = c('ID', paste(indiv, samp, colnames(expression)[-1], sep='_'))
      if(is.null(pseudo_bulk_tcells_allGenes)) pseudo_bulk_tcells_allGenes = expression
      else pseudo_bulk_tcells_allGenes = merge(pseudo_bulk_tcells_allGenes, expression, by='ID')
    }
  }
}
rownames(pseudo_bulk_tcells_allGenes) = pseudo_bulk_tcells_allGenes$ID
pseudo_bulk_tcells_allGenes = pseudo_bulk_tcells_allGenes[,-1]

pseudo_bulk_tcells = pseudo_bulk_tcells_allGenes[metabolic_genes[metabolic_genes%in%rownames(pseudo_bulk_tcells_allGenes)],]
for(col_name in colnames(pseudo_bulk_tcells)){
  pseudo_bulk_tcells[,col_name] = as.numeric(pseudo_bulk_tcells[,col_name])
}



# -----------------------
# --- Reaction fluxes ---
# -----------------------

# Read FBA data:
fba_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/normal_FBA.csv', row.names=1, check.names = FALSE)
# Sigmoid Transform fluxes:
fba_fluxes = 2 * ((1 / (1 + exp(-fba_fluxes_raw))) - 0.5)
fba_fluxes_abs = abs(fba_fluxes) # as gene expression is just abundance and not direction
# Create metadata:
split_colnames = strsplit(colnames(fba_fluxes), '_')
samples = c()
cell_types = c()
individuals = c()
for(row_split in split_colnames){
  individuals = c(individuals, row_split[1])
  samples = c(samples, row_split[2])
  cell_types = c(cell_types, row_split[3])
}
metadata = data.frame(individual=individuals, sample=samples, cell_type=cell_types,
                      row.names=colnames(fba_fluxes))
more_meta = read.csv('./0Data/scRNAseq/CRC_atlas/expression_data/metadata.csv', row.names=1)
states = c()
for(samp in metadata$sample){
  states = c(states, more_meta[samp, 'Sample.Source'])
}
metadata$state = states



# --------------------------------------------------------------
# --- Reaction fluxes and presence, only reactions with GPRs ---
# --------------------------------------------------------------

# 1. GPRs:
gprs = read.csv('./GENERAL/utility_data/HumanGEM-1.8.0_consistent_GPRs.txt', sep='\t', header=FALSE,
                row.names = 1)
# 1.1. Get reactions with GPRs
rxns_with_gpr = rownames(gprs)[gprs$V2 != ''] # 7 348
# 1.2. From the reactions with GPRs, get only those with at least one gene in the rna dataset
# 1.2.1. Map rna's gene symbols to ensemble ids
pseudo_bulk_tcells_ensemble = pseudo_bulk_tcells
rownames(pseudo_bulk_tcells_ensemble) = humangem_genemapping[rownames(pseudo_bulk_tcells_ensemble)]
pseudo_bulk_tcells_ensemble = t(scale(t(pseudo_bulk_tcells_ensemble)))
# 1.2.2. Check reactions with genes present in rna dataset:
rxns_keep = c()
for(rxn in rxns_with_gpr){
  rxn_genes = grep('ENSG', strsplit(gprs[rxn, 'V2'], ' ')[[1]], value=T)
  present_rxn_genes = rxn_genes %in% rownames(pseudo_bulk_tcells_ensemble)
  if(sum(present_rxn_genes)>0) rxns_keep = c(rxns_keep, rxn)
}
# 1.2.2. 6 297 reactions will be kept

# 2. Get fluxes dataset only with these reactions
fba_fluxes_gprs = fba_fluxes[rxns_keep, ]
fba_fluxes_abs_gprs = fba_fluxes_abs[rxns_keep, ]

# 3. Get reaction activity with all reactions
fba_activity = fba_fluxes
fba_activity[fba_activity>0] = 1
fba_activity[fba_activity<0] = -1

# 4. Get reaction activity only with these reactions
fba_activity_gprs = fba_fluxes_gprs
fba_activity_gprs[fba_activity_gprs>0] = 1
fba_activity_gprs[fba_activity_gprs<0] = -1

# 5. Get reaction presence only with the reactions
rxn_presence = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/general/reaction_presence_Tcells.csv',
                        row.names=1, header=TRUE, check.names=FALSE)
rxn_presence_gprs = rxn_presence[rxns_keep, ]



# --------------------------------------
# --- Save environment for rmarkdown ---
# --------------------------------------

save.image('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/general/reports/NormalMatched_env.RData')



# -----------------------------------------------
# --- UMAP of scRNAseq dataset with all genes ---
# -----------------------------------------------

# UMAP of all cells:
scRNAseq_allCells_allGenes = Seurat::DimPlot(tcells, group.by='models_annotation', label=FALSE, pt.size=.5) +
  ggplot2::scale_colour_manual(values=ct_colors) + ggplot2::ggtitle('scRNAseq | All Genes') +
  ggplot2::theme(legend.position = 'bottom')



# --------------------------------------------------------------
# --- UMAP of scRNAseq dataset with only the metabolic genes ---
# --------------------------------------------------------------

# PCA on all cells:
tcells_metabolic = Seurat::NormalizeData(tcells_metabolic, assay='RNA')
tcells_metabolic = Seurat::ScaleData(tcells_metabolic)
tcells_metabolic = Seurat::RunPCA(tcells_metabolic, assay='integrated')
tcells_metabolic = Seurat::RunUMAP(tcells_metabolic, dims = 1:30, verbose = FALSE)
invisible(gc())
scRNAseq_allCells_metabolicGenes = Seurat::DimPlot(tcells_metabolic, group.by='models_annotation', label=FALSE, pt.size=.5) +
  ggplot2::scale_colour_manual(values=ct_colors) + ggplot2::ggtitle('scRNAseq | Metabolic Genes') +
  ggplot2::theme(legend.position = 'none')





# ------------------------------------------------------------------
# --- UMAP of pseudo-bulk dataset with only the metabolic genes ---
# -----------------------------------------------------------------

pseudo_bulk_tcells_seurat = SeuratObject::CreateSeuratObject(pseudo_bulk_tcells,
                                                             meta.data = metadata[colnames(pseudo_bulk_tcells),])

pseudo_bulk_tcells_seurat = Seurat::NormalizeData(pseudo_bulk_tcells_seurat)
pseudo_bulk_tcells_seurat = Seurat::ScaleData(pseudo_bulk_tcells_seurat)
pseudo_bulk_tcells_seurat = Seurat::RunPCA(pseudo_bulk_tcells_seurat, assay='RNA',
                                           features=rownames(pseudo_bulk_tcells))
pseudo_bulk_tcells_seurat = Seurat::RunUMAP(pseudo_bulk_tcells_seurat, dims = 1:30, verbose = FALSE)
invisible(gc())
pseudoBulk_allCellTypes_metabolicGenes = Seurat::DimPlot(pseudo_bulk_tcells_seurat, group.by='cell_type', label=FALSE, pt.size=1.5) +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::ggtitle('pseudo-bulk | Metabolic Genes') + ggplot2::theme(legend.position='bottom')





# -------------------------------------------------------------------------
# --- UMAP of reaction fluxes dataset with only the reactions with GPRs ---
# -------------------------------------------------------------------------

fba_fluxes_gprs_tcells_seurat = SeuratObject::CreateSeuratObject(fba_fluxes_gprs,
                                                                 meta.data=metadata[colnames(fba_fluxes_gprs),])

fba_fluxes_gprs_tcells_seurat = Seurat::ScaleData(fba_fluxes_gprs_tcells_seurat)
fba_fluxes_gprs_tcells_seurat = Seurat::RunPCA(fba_fluxes_gprs_tcells_seurat, assay='RNA',
                                               features=rownames(fba_fluxes_gprs))
fba_fluxes_gprs_tcells_seurat = Seurat::RunUMAP(fba_fluxes_gprs_tcells_seurat, dims = 1:30, verbose = FALSE)
invisible(gc())
fluxes_allCellTypes_GPRrxns = Seurat::DimPlot(fba_fluxes_gprs_tcells_seurat, group.by='cell_type', label=FALSE, pt.size=1.5) +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::ggtitle('reaction fluxes | reactions with GPRs') + ggplot2::theme(legend.position='bottom')



# ----------------------------------------------------------------------------------
# --- UMAP of reaction absolute fluxes dataset with only the reactions with GPRs ---
# ----------------------------------------------------------------------------------

fba_fluxes_abs_gprs_tcells_seurat = SeuratObject::CreateSeuratObject(fba_fluxes_abs_gprs,
                                                                     meta.data=metadata[colnames(fba_fluxes_abs_gprs),])

fba_fluxes_abs_gprs_tcells_seurat = Seurat::ScaleData(fba_fluxes_abs_gprs_tcells_seurat)
fba_fluxes_abs_gprs_tcells_seurat = Seurat::RunPCA(fba_fluxes_abs_gprs_tcells_seurat, assay='RNA',
                                               features=rownames(fba_fluxes_abs_gprs))
fba_fluxes_abs_gprs_tcells_seurat = Seurat::RunUMAP(fba_fluxes_abs_gprs_tcells_seurat, dims = 1:30, verbose = FALSE)
invisible(gc())
fluxesAbs_allCellTypes_GPRrxns = Seurat::DimPlot(fba_fluxes_abs_gprs_tcells_seurat, group.by='cell_type', label=FALSE, pt.size=1.5) +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::ggtitle('absolute reaction fluxes | reactions with GPRs') + ggplot2::theme(legend.position='bottom')



# ---------------------------------------------------------------------------
# --- UMAP of reaction presence dataset with only the reactions with GPRs ---
# ---------------------------------------------------------------------------

rxn_presence_gprs_tcells_seurat = SeuratObject::CreateSeuratObject(rxn_presence_gprs,
                                                                     meta.data=metadata[colnames(rxn_presence_gprs),])

rxn_presence_gprs_tcells_seurat = Seurat::ScaleData(rxn_presence_gprs_tcells_seurat)
rxn_presence_gprs_tcells_seurat = Seurat::RunPCA(rxn_presence_gprs_tcells_seurat, assay='RNA',
                                                   features=rownames(rxn_presence_gprs))
rxn_presence_gprs_tcells_seurat = Seurat::RunUMAP(rxn_presence_gprs_tcells_seurat, dims = 1:30, verbose = FALSE)
invisible(gc())
rxn_presence_allCellTypes_GPRrxns = Seurat::DimPlot(rxn_presence_gprs_tcells_seurat, group.by='cell_type', label=FALSE, pt.size=1.5) +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::ggtitle('reaction presence | reactions with GPRs') + ggplot2::theme(legend.position='bottom')



# -------------------------------
# --- Join the UMAPs together ---
# -------------------------------

(scRNAseq_allCells_allGenes + ggplot2::theme(legend.position='none') | 
    scRNAseq_allCells_metabolicGenes + ggplot2::theme(legend.position='none')) / 
  (pseudoBulk_allCellTypes_metabolicGenes + ggplot2::theme(legend.position='none') | 
     rxn_presence_allCellTypes_GPRrxns + ggplot2::theme(legend.position='none')) /
  (fluxes_allCellTypes_GPRrxns | fluxesAbs_allCellTypes_GPRrxns + ggplot2::theme(legend.position='none'))

(scRNAseq_allCells_allGenes + ggplot2::theme(legend.position='none') + ggplot2::ggtitle('A') | 
    scRNAseq_allCells_metabolicGenes + ggplot2::theme(legend.position='none') + ggplot2::ggtitle('B')) / 
  (pseudoBulk_allCellTypes_metabolicGenes + ggplot2::theme(legend.position='none') + ggplot2::ggtitle('C') | 
     rxn_presence_allCellTypes_GPRrxns + ggplot2::theme(legend.position='none') + ggplot2::ggtitle('D')) /
  (fluxes_allCellTypes_GPRrxns + ggplot2::ggtitle('E') | fluxesAbs_allCellTypes_GPRrxns + ggplot2::theme(legend.position='none') + ggplot2::ggtitle('F'))





# ------------------------
# --- Machine Learning ---
# ------------------------


# 0. Train and test samples:
train_samples = sample(colnames(pseudo_bulk_tcells), 120)
test_samples = colnames(pseudo_bulk_tcells)[!colnames(pseudo_bulk_tcells)%in%train_samples]
# Distribution of cell-types in train samples
table(metadata[train_samples, 'cell_type']) / table(metadata[, 'cell_type']) * 100
jsonlite::write_json(list(train=train_samples, test=test_samples),
                     './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/train_test_samples.json')
# Where data will be stored:
ml_res = data.frame(MCC=rep(NA, 10), row.names=c('pseudoBulk_allGenes', 'pseudoBulk_metabolicGenes',
                                                 'fbaFluxes_allRxns', 'fbaFluxes_GPRs',
                                                 'fbaAbsFluxes_allRxns', 'fbaAbsFluxes_GPRs',
                                                 'fbaActivity_allRxns', 'fbaActivity_GPRs',
                                                 'rxnPresence_allRxns', 'rxnPresence_GPRs'))
ml_varImp = list()


# 1. Pseudo-bulk gene expression (all genes)
# 1.1. Dataset for machine learning:
ml_pseudoAll_data = t(pseudo_bulk_tcells_allGenes)
saveRDS(ml_pseudoAll_data, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/1_pseudoBulkAllGenes_data.RData')

# 1.2. Feature Selection
# 1.2.1. Remove low variance reactions
genes_to_remove = caret::nearZeroVar(ml_pseudoAll_data[train_samples,])

# 1.3. Random Forest
# 1.3.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_pseudoAll = caret::train(x = ml_pseudoAll_data[train_samples,-genes_to_remove], 
                                     y = as.factor(metadata[train_samples, 'cell_type']),
                                     method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_pseudoAll, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/1_pseudoBulkAllGenes_model.RData')
# 1.3.2. Predict test data:
predictions_rf_pseudoAll = predict(cv_model_rf_pseudoAll, ml_pseudoAll_data[test_samples,])
# 1.3.3. Measure prediction capacity
rf_pseudoAll_cm = caret::confusionMatrix(predictions_rf_pseudoAll, as.factor(metadata[test_samples, 'cell_type']))$table
write.csv(rf_pseudoAll_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/1_pseudoBulkAllGenes_confusionMatrix.csv')
rf_pseudoAll_mcc = yardstick::mcc(data.frame(pred=predictions_rf_pseudoAll,
                                             obs=as.factor(metadata[test_samples, 'cell_type'])),
                                  'obs', 'pred')
ml_res['pseudoBulk_allGenes', 'MCC'] = rf_pseudoAll_mcc$.estimate
# 1.3.4. Variable Importance:
ml_varImp[['pseudoBulk_allGenes']] = caret::varImp(cv_model_rf_pseudoAll)$importance


# 2. Pseudo-bulk gene expression (metabolic genes)
# 2.1. Dataset for machine learning:
ml_pseudo_data = t(pseudo_bulk_tcells)
saveRDS(ml_pseudo_data, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/2_pseudoBulkMetabolicGenes_data.RData')

# 2.2. Feature Selection
# 2.2.1. Remove low variance reactions
genes_to_remove = caret::nearZeroVar(ml_pseudo_data[train_samples,])

# 2.3. Random Forest
# 2.3.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_pseudo = caret::train(x = ml_pseudo_data[train_samples,-genes_to_remove], 
                                  y = as.factor(metadata[train_samples, 'cell_type']),
                                  method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_pseudo, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/2_pseudoBulkMetabolicGenes_model.RData')
# 2.3.2. Predict test data:
predictions_rf_pseudo = predict(cv_model_rf_pseudo, ml_pseudo_data[test_samples,])
# 2.3.3. Measure prediction capacity
rf_pseudo_cm = caret::confusionMatrix(predictions_rf_pseudo, as.factor(metadata[test_samples, 'cell_type']))$table
write.csv(rf_pseudo_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/2_pseudoBulkMetabolicGenes_confusionMatrix.csv')
rf_pseudo_mcc = yardstick::mcc(data.frame(pred=predictions_rf_pseudo,
                                           obs=as.factor(metadata[test_samples, 'cell_type'])),
                                'obs', 'pred')
ml_res['pseudoBulk_metabolicGenes', 'MCC'] = rf_pseudo_mcc$.estimate
# 2.3.4. Variable Importance:
ml_varImp[['pseudoBulk_metabolicGenes']] = caret::varImp(cv_model_rf_pseudo)$importance


# 3. Reaction fluxes (all reactions)
# 3.1. Dataset for machine learning:
ml_fbaFluxes_data = t(fba_fluxes[, colnames(pseudo_bulk_tcells)])
saveRDS(ml_fbaFluxes_data, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/3_fbaFluxesAllRxns_data.RData')

# 3.2. Feature Selection
# 3.2.1. Remove low variance reactions
rxns_to_remove = caret::nearZeroVar(ml_fbaFluxes_data[train_samples,])

# 3.3. Random Forest
# 3.3.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_fbaFluxes = caret::train(x = ml_fbaFluxes_data[train_samples,-rxns_to_remove], 
                                     y = as.factor(metadata[train_samples, 'cell_type']),
                                     method = 'rf',
                                     trControl = fitControl)
saveRDS(cv_model_rf_fbaFluxes, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/3_fbaFluxesAllRxns_model.RData')
# 3.3.2. Predict test data:
predictions_rf_fbaFluxes = predict(cv_model_rf_fbaFluxes, ml_fbaFluxes_data[test_samples,])
# 3.3.3. Measure prediction capacity
rf_fbaFluxes_cm = caret::confusionMatrix(predictions_rf_fbaFluxes, as.factor(metadata[test_samples, 'cell_type']))$table
write.csv(rf_fbaFluxes_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/3_fbaFluxesAllRxns_confusionMatrix.csv')
rf_fbaFluxes_mcc = yardstick::mcc(data.frame(pred=predictions_rf_fbaFluxes,
                                             obs=as.factor(metadata[test_samples, 'cell_type'])),
                                  'obs', 'pred')
ml_res['fbaFluxes_allRxns', 'MCC'] = rf_fbaFluxes_mcc$.estimate
# 3.3.4. Variable Importance:
ml_varImp[['fbaFluxes_allRxns']] = caret::varImp(cv_model_rf_fbaFluxes)$importance


# 4. Reaction fluxes (GPRs reactions)
# 4.1. Dataset for machine learning:
ml_fbaFluxesGPRs_data = t(fba_fluxes_gprs[, colnames(pseudo_bulk_tcells)])
saveRDS(ml_fbaFluxesGPRs_data, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/4_fbaFluxesGPRs_data.RData')

# 4.2. Feature Selection
# 4.2.1. Remove low variance reactions
rxns_to_remove = caret::nearZeroVar(ml_fbaFluxesGPRs_data[train_samples,])

# 4.3. Random Forest
# 4.3.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_fbaFluxesGPRs = caret::train(x = ml_fbaFluxesGPRs_data[train_samples,-rxns_to_remove],
                                         y = as.factor(metadata[train_samples, 'cell_type']),
                                         method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_fbaFluxesGPRs, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/4_fbaFluxesGPRs_model.RData')
# 4.3.2. Predict test data:
predictions_rf_fbaFluxesGPRs = predict(cv_model_rf_fbaFluxesGPRs, ml_fbaFluxesGPRs_data[test_samples,])
# 4.3.3. Measure prediction capacity
rf_fbaFluxesGPRs_cm = caret::confusionMatrix(predictions_rf_fbaFluxesGPRs, as.factor(metadata[test_samples, 'cell_type']))$table
write.csv(rf_fbaFluxesGPRs_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/4_fbaFluxesGPRs_confusionMatrix.csv')
rf_fbaFluxesGPRs_mcc = yardstick::mcc(data.frame(pred=predictions_rf_fbaFluxesGPRs,
                                                 obs=as.factor(metadata[test_samples, 'cell_type'])),
                                      'obs', 'pred')
ml_res['fbaFluxes_GPRs', 'MCC'] = rf_fbaFluxesGPRs_mcc$.estimate
# 4.3.4. Variable Importance:
ml_varImp[['fbaFluxes_GPRs']] = caret::varImp(cv_model_rf_fbaFluxesGPRs)$importance


# 5. Reaction absolute fluxes (all reactions)
# 5.1. Dataset for machine learning:
ml_fbaFluxesAbs_data = t(fba_fluxes_abs[, colnames(pseudo_bulk_tcells)])
saveRDS(ml_fbaFluxesAbs_data, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/5_fbaAbsFluxesAllRxns_data.RData')

# 5.2. Feature Selection
# 5.2.1. Remove low variance reactions
rxns_to_remove = caret::nearZeroVar(ml_fbaFluxesAbs_data[train_samples,])

# 5.3. Random Forest
# 5.3.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_fbaFluxesAbs = caret::train(x = ml_fbaFluxesAbs_data[train_samples,-rxns_to_remove],
                                         y = as.factor(metadata[train_samples, 'cell_type']),
                                         method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_fbaFluxesAbs, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/5_fbaAbsFluxesAllRxns_model.RData')
# 5.3.2. Predict test data:
predictions_rf_fbaFluxesAbs = predict(cv_model_rf_fbaFluxesAbs, ml_fbaFluxesAbs_data[test_samples,])
# 5.3.3. Measure prediction capacity
rf_fbaFluxesAbs_cm = caret::confusionMatrix(predictions_rf_fbaFluxesAbs, as.factor(metadata[test_samples, 'cell_type']))$table
write.csv(rf_fbaFluxesAbs_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/5_fbaAbsFluxesAllRxns_confusionMatrix.csv')
rf_fbaFluxesAbs_mcc = yardstick::mcc(data.frame(pred=predictions_rf_fbaFluxesAbs,
                                                 obs=as.factor(metadata[test_samples, 'cell_type'])),
                                      'obs', 'pred')
ml_res['fbaAbsFluxes_allRxns', 'MCC'] = rf_fbaFluxesAbs_mcc$.estimate
# 4.3.4. Variable Importance:
ml_varImp[['fbaAbsFluxes_allRxns']] = caret::varImp(cv_model_rf_fbaFluxesAbs)$importance


# 6. Reaction absolute fluxes (GPRs reactions)
# 6.1. Dataset for machine learning:
ml_fbaFluxesAbsGPRs_data = t(fba_fluxes_abs_gprs[, colnames(pseudo_bulk_tcells)])
saveRDS(ml_fbaFluxesAbsGPRs_data, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/6_fbaAbsFluxesGPRs_data.RData')

# 6.2. Feature Selection
# 6.2.1. Remove low variance reactions
rxns_to_remove = caret::nearZeroVar(ml_fbaFluxesAbsGPRs_data[train_samples,])

# 6.3. Random Forest
# 6.3.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_fbaFluxesAbsGPRs = caret::train(x = ml_fbaFluxesAbsGPRs_data[train_samples,-rxns_to_remove],
                                            y = as.factor(metadata[train_samples, 'cell_type']),
                                            method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_fbaFluxesAbsGPRs, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/6_fbaAbsFluxesGPRs_model.RData')
# 6.3.2. Predict test data:
predictions_rf_fbaFluxesAbsGPRs = predict(cv_model_rf_fbaFluxesAbsGPRs, ml_fbaFluxesAbsGPRs_data[test_samples,])
# 6.3.3. Measure prediction capacity
rf_fbaFluxesAbsGPRs_cm = caret::confusionMatrix(predictions_rf_fbaFluxesAbsGPRs, as.factor(metadata[test_samples, 'cell_type']))$table
write.csv(rf_fbaFluxesAbsGPRs_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/6_fbaAbsFluxesGPRs_confusionMatrix.csv')
rf_fbaFluxesAbsGPRs_mcc = yardstick::mcc(data.frame(pred=predictions_rf_fbaFluxesAbsGPRs,
                                                    obs=as.factor(metadata[test_samples, 'cell_type'])),
                                         'obs', 'pred')
ml_res['fbaAbsFluxes_GPRs', 'MCC'] = rf_fbaFluxesAbsGPRs_mcc$.estimate
# 4.3.4. Variable Importance:
ml_varImp[['fbaAbsFluxes_GPRs']] = caret::varImp(cv_model_rf_fbaFluxesAbsGPRs)$importance


# 7. Reaction Activity (All reactions)
# 7.1. Dataset for machine learning:
ml_fbaActivity_data = t(fba_activity[, colnames(pseudo_bulk_tcells)])
saveRDS(ml_fbaActivity_data, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/7_fbaActivityAllRxns_data.RData')

# 7.2. Feature Selection
# 7.2.1. Remove low variance reactions
rxns_to_remove = caret::nearZeroVar(ml_fbaActivity_data[train_samples,])

# 7.3. Random Forest
# 7.3.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_fbaActivity = caret::train(x = ml_fbaActivity_data[train_samples,-rxns_to_remove],
                                            y = as.factor(metadata[train_samples, 'cell_type']),
                                            method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_fbaActivity, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/7_fbaActivityAllRxns_model.RData')
# 7.3.2. Predict test data:
predictions_rf_fbaActivity = predict(cv_model_rf_fbaActivity, ml_fbaActivity_data[test_samples,])
# 7.3.3. Measure prediction capacity
rf_fbaActivity_cm = caret::confusionMatrix(predictions_rf_fbaActivity, as.factor(metadata[test_samples, 'cell_type']))$table
write.csv(rf_fbaActivity_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/7_fbaActivityAllRxns_confusionMatrix.csv')
rf_fbaActivity_mcc = yardstick::mcc(data.frame(pred=predictions_rf_fbaActivity,
                                                    obs=as.factor(metadata[test_samples, 'cell_type'])),
                                         'obs', 'pred')
ml_res['fbaActivity_allRxns', 'MCC'] = rf_fbaActivity_mcc$.estimate
# 4.3.4. Variable Importance:
ml_varImp[['fbaActivity_allRxns']] = caret::varImp(cv_model_rf_fbaActivity)$importance


# 8. Reaction Activity (GPRs reactions)
# 8.1. Dataset for machine learning:
ml_fbaActivityGPRs_data = t(fba_activity_gprs[, colnames(pseudo_bulk_tcells)])
saveRDS(ml_fbaActivityGPRs_data, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/8_fbaActivityGPRs_data.RData')

# 8.2. Feature Selection
# 8.2.1. Remove low variance reactions
rxns_to_remove = caret::nearZeroVar(ml_fbaActivityGPRs_data[train_samples,])

# 8.3. Random Forest
# 8.3.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_fbaActivityGPRs = caret::train(x = ml_fbaActivityGPRs_data[train_samples,-rxns_to_remove],
                                           y = as.factor(metadata[train_samples, 'cell_type']),
                                           method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_fbaActivityGPRs, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/8_fbaActivityGPRs_model.RData')
# 8.3.2. Predict test data:
predictions_rf_fbaActivityGPRs = predict(cv_model_rf_fbaActivityGPRs, ml_fbaActivityGPRs_data[test_samples,])
# 8.3.3. Measure prediction capacity
rf_fbaActivityGPRs_cm = caret::confusionMatrix(predictions_rf_fbaActivityGPRs, as.factor(metadata[test_samples, 'cell_type']))$table
write.csv(rf_fbaActivityGPRs_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/8_fbaActivityGPRs_confusionMatrix.csv')
rf_fbaActivityGPRs_mcc = yardstick::mcc(data.frame(pred=predictions_rf_fbaActivityGPRs,
                                                   obs=as.factor(metadata[test_samples, 'cell_type'])),
                                        'obs', 'pred')
ml_res['fbaActivity_GPRs', 'MCC'] = rf_fbaActivityGPRs_mcc$.estimate
# 4.3.4. Variable Importance:
ml_varImp[['fbaActivity_GPRs']] = caret::varImp(cv_model_rf_fbaActivityGPRs)$importance


# 9. Reaction Presence (All reactions)
# 9.1. Dataset for machine learning:
ml_rxnPresence_data = t(rxn_presence[, colnames(pseudo_bulk_tcells)])
saveRDS(ml_rxnPresence_data, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/9_rxnPresenceAllRxns_data.RData')

# 9.2. Feature Selection
# 9.2.1. Remove low variance reactions
rxns_to_remove = caret::nearZeroVar(ml_rxnPresence_data[train_samples,])

# 9.3. Random Forest
# 9.3.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_rxnPresence = caret::train(x = ml_rxnPresence_data[train_samples,-rxns_to_remove],
                                           y = as.factor(metadata[train_samples, 'cell_type']),
                                           method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_rxnPresence, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/9_rxnPresenceAllRxns_model.RData')
# 9.3.2. Predict test data:
predictions_rf_rxnPresence = predict(cv_model_rf_rxnPresence, ml_rxnPresence_data[test_samples,])
# 9.3.3. Measure prediction capacity
rf_rxnPresence_cm = caret::confusionMatrix(predictions_rf_rxnPresence, as.factor(metadata[test_samples, 'cell_type']))$table
write.csv(rf_rxnPresence_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/9_rxnPresenceAllRxns_confusionMatrix.csv')
rf_rxnPresence_mcc = yardstick::mcc(data.frame(pred=predictions_rf_rxnPresence,
                                               obs=as.factor(metadata[test_samples, 'cell_type'])),
                                    'obs', 'pred')
ml_res['rxnPresence_allRxns', 'MCC'] = rf_rxnPresence_mcc$.estimate
# 4.3.4. Variable Importance:
ml_varImp[['rxnPresence_allRxns']] = caret::varImp(cv_model_rf_rxnPresence)$importance


# 10. Reaction Presence (All reactions)
# 10.1. Dataset for machine learning:
ml_rxnPresenceGPRs_data = t(rxn_presence_gprs[, colnames(pseudo_bulk_tcells)])
saveRDS(ml_rxnPresenceGPRs_data, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/10_rxnPresenceGPRs_data.RData')

# 10.2. Feature Selection
# 10.2.1. Remove low variance reactions
rxns_to_remove = caret::nearZeroVar(ml_rxnPresenceGPRs_data[train_samples,])

# 10.3. Random Forest
# 10.3.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_rxnPresenceGPRs = caret::train(x = ml_rxnPresenceGPRs_data[train_samples,-rxns_to_remove],
                                           y = as.factor(metadata[train_samples, 'cell_type']),
                                           method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_rxnPresenceGPRs, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/10_rxnPresenceGPRs_model.RData')
# 10.3.2. Predict test data:
predictions_rf_rxnPresenceGPRs = predict(cv_model_rf_rxnPresenceGPRs, ml_rxnPresenceGPRs_data[test_samples,])
# 10.3.3. Measure prediction capacity
rf_rxnPresenceGPRs_cm = caret::confusionMatrix(predictions_rf_rxnPresenceGPRs, as.factor(metadata[test_samples, 'cell_type']))$table
write.csv(rf_rxnPresenceGPRs_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/10_rxnPresenceGPRs_confusionMatrix.csv')
rf_rxnPresenceGPRs_mcc = yardstick::mcc(data.frame(pred=predictions_rf_rxnPresenceGPRs,
                                               obs=as.factor(metadata[test_samples, 'cell_type'])),
                                    'obs', 'pred')
ml_res['rxnPresence_GPRs', 'MCC'] = rf_rxnPresenceGPRs_mcc$.estimate
# 4.3.4. Variable Importance:
ml_varImp[['rxnPresence_GPRs']] = caret::varImp(cv_model_rf_rxnPresenceGPRs)$importance


# 11. Store MCC metric values and variable importance list:
write.csv(ml_res, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/MCC.csv')
jsonlite::write_json(ml_varImp, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/varImp.json')





# ---------------------------------------------------------------------
# --- Over-Representation Analysis (ORA) on the variable importance ---
# ---------------------------------------------------------------------

# 0.1. Load annotations and all variable importances:
HumanGEM_subsystems_rxns = jsonlite::read_json('./GENERAL/utility_data/subsystems_reactions_mapping.json',
                                               simplifyVector = TRUE)
for(path in names(HumanGEM_subsystems_rxns)){
  HumanGEM_subsystems_rxns[[path]] = HumanGEM_subsystems_rxns[[path]][HumanGEM_subsystems_rxns[[path]]%in%rownames(rxn_presence_gprs)]
}

paths_annots = c()
for(path in names(HumanGEM_subsystems_rxns)){
  path_rxns = cbind(rep(path, length(HumanGEM_subsystems_rxns[[path]])),
                    HumanGEM_subsystems_rxns[[path]])
  if(length(paths_annots)==0) paths_annots = path_rxns
  else paths_annots = rbind(paths_annots, path_rxns)
}

var_imp = jsonlite::read_json('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ml/varImp.json',
                              simplifyVector = TRUE)

# 1. Reactions Presence Data:
# 1.1. Get ranked list of reactions:
varImp_rxnpresence = as.numeric(var_imp$rxnPresence_GPRs$Overall)
names(varImp_rxnpresence) = rownames(var_imp$rxnPresence_GPRs)
varImp_rxnpresence = sort(varImp_rxnpresence, decreasing = TRUE)
# 1.2. Perform the enrichment analysis:
gsea_rxnPresence = clusterProfiler::GSEA(varImp_rxnpresence, TERM2GENE=paths_annots,
                                         pAdjustMethod='fdr', pvalueCutoff=1)
# 1.3. Heatmap of enrichment analysis' pathways
bile_acid_varImp = names(sort(varImp_rxnpresence[HumanGEM_subsystems_rxns$`Bile acid biosynthesis`], decreasing=TRUE)[sort(varImp_rxnpresence[HumanGEM_subsystems_rxns$`Bile acid biosynthesis`], decreasing=TRUE) > 6])
rxnPresence_cols = c('lightgrey', 'black')
names(rxnPresence_cols) = c('0', '1')
ComplexHeatmap::Heatmap(as.matrix(rxn_presence_gprs[bile_acid_varImp,]),
                        row_order = 1:12,
                        column_split = metadata[colnames(rxn_presence_gprs),]$cell_type,
                        show_column_names = FALSE,
                        col = rxnPresence_cols,
                        top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata[colnames(rxn_presence_gprs),]$cell_type,
                                                                         state=metadata[colnames(rxn_presence_gprs),]$state,
                                                                         col=list(cell_type=ct_colors, state=state_colors)))




