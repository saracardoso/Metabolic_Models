ct_colors = c('#b30000', '#999999', '#cc79a7', '#b3c9e6', '#4e84c4', '#c3d7a4', '#52854c', '#fdb981',
              '#d16103', '#293352')
names(ct_colors) = c('Cytotoxic CD8 Tcells', 'Follicular CD4 Tcells', 'IL17+ CD4 Tcells',
                     'Memory CD4 Tcells', 'Memory CD8 Tcells', 'Naive CD4 Tcells', 'Naive CD8 Tcells',
                     'Proliferative CD4 Tcells', 'Proliferative CD8 Tcells',  'Regulatory CD4 Tcells')
state_colors = c('#52854c', '#b30000')
names(state_colors) = c('Normal Matched', 'Tumour')

# ---------------------------
# --- Metadata and others ---
# ---------------------------

metadata = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/metadata.csv', row.names=1)
genes_pathways = read.csv('./GENERAL/utility_data/genes_subsystems_mapping.csv')


# ----------------------------------------------
# --- scRNAseq data used for models analysed ---
# ----------------------------------------------

# Read scRNAseq CRC atlas:
tcells = SeuratDisk::LoadH5Seurat('/home/scardoso/Documents/PhD/CRC_ATLAS/2_annotation/results_Tcells/datasets/Tcells_finalAnnots.h5Seurat')

# Get only individuals that were used in the models:
indivs_keep = rownames(tcells[[]])[tcells$patient%in%c('31', '32', '33', '35', 'KUL01', 'KUL19',
                                                       'KUL21', 'SMC01', 'SMC04', 'SMC06', 'SMC07',
                                                       'SMC08', 'SMC10')]
tcells = subset(tcells, cells=indivs_keep)
invisible(gc())

# Get cell-type annotation equal to the models:
tcells[['models_annotation']] = tcells$Annotation_Level_4

# Get only the cells from the models that ended up being analysed:
cells_keep = c()
for(samp in unique(metadata$sample)){
  ct = metadata$cell_type[metadata$sample==samp]
  cells_keep = c(cells_keep,
                 rownames(tcells[[]])[tcells$sample==samp & tcells$models_annotation%in%ct])
}
tcells = subset(tcells, cells=cells_keep)
invisible(gc())
# Store the scRNAseq data:
SeuratDisk::SaveH5Seurat(tcells,
                         './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/tcells_allGenes.h5Seurat')



# --------------------------------------------------------------------
# --- scRNAseq data used for models analysed, only metabolic genes ---
# --------------------------------------------------------------------

metabolic_genes = unique(genes_pathways$gene_symbol)
metabolic_genes_scrnaseq = metabolic_genes[metabolic_genes%in%rownames(tcells@assays$RNA@counts)]
tcells_metabolic = subset(tcells, features = metabolic_genes_scrnaseq)
invisible(gc())
# Store the scRNAseq data:
SeuratDisk::SaveH5Seurat(tcells_metabolic,
                         './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/tcells_metabolicGenes.h5Seurat')





# ---------------------------------------------------------------
# --- pseudo-bulk data (CPMs, GASs, RASs) for models analysed ---
# ---------------------------------------------------------------

data_folder = './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched'
cpm_counts_tcells = c()
tas_counts_tcells = c()
ras_counts_tcells = c()
for(samp in unique(metadata$sample)){
  message(samp)
  indiv = unique(metadata$individual[metadata$sample==samp])
  cpm_counts = read.csv(paste(data_folder, '/', indiv, '/data_info/', samp, '_CPM.csv', sep=''),
                        row.names=1, check.names = FALSE)
  tas_counts = read.csv(paste(data_folder, '/', indiv, '/data_info/', samp, '_TAS.csv', sep=''),
                        row.names=1, check.names = FALSE)
  ras_counts = read.csv(paste(data_folder, '/', indiv, '/data_info/', samp, '_RAS.csv', sep=''),
                        row.names=1, check.names = FALSE)
  
  colnames(cpm_counts) = paste(indiv, samp, colnames(cpm_counts), sep='_')
  colnames(tas_counts) = paste(indiv, samp, colnames(tas_counts), sep='_')
  colnames(ras_counts) = paste(indiv, samp, colnames(ras_counts), sep='_')
  
  cpm_counts = cpm_counts[, colnames(cpm_counts)%in%rownames(metadata)]
  tas_counts = tas_counts[, colnames(tas_counts)%in%rownames(metadata)]
  ras_counts = ras_counts[, colnames(ras_counts)%in%rownames(metadata)]
  
  if(is.null(cpm_counts_tcells)){
    cpm_counts_tcells = cpm_counts
    tas_counts_tcells = tas_counts
    ras_counts_tcells = ras_counts
  }
  else{
    cpm_counts_tcells = merge(cpm_counts_tcells, cpm_counts, by="row.names", all.x=TRUE)
    rownames(cpm_counts_tcells) = cpm_counts_tcells[,1]
    cpm_counts_tcells = cpm_counts_tcells[,-1]
    tas_counts_tcells = cbind(tas_counts_tcells, tas_counts)
    ras_counts_tcells = cbind(ras_counts_tcells, ras_counts)
  }
}
cpm_counts_tcells[is.na(cpm_counts_tcells)] = 0
tas_counts_tcells[is.na(tas_counts_tcells)] = 0
ras_counts_tcells[is.na(ras_counts_tcells)] = 0

genes_map = unique(genes_pathways[,c('gene', 'gene_symbol')])
genes_with_info_ensemble = genes_map[genes_map$gene_symbol%in%rownames(cpm_counts_tcells), ]

cpm_counts_tcells = cpm_counts_tcells[genes_with_info_ensemble$gene_symbol,]
tas_counts_tcells = tas_counts_tcells[genes_with_info_ensemble$gene, ]
rownames(tas_counts_tcells) = genes_with_info_ensemble$gene_symbol

write.csv(cpm_counts_tcells,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/CPMs.csv')
write.csv(tas_counts_tcells,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/GASs.csv')
write.csv(ras_counts_tcells,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/RASs.csv')




# ------------------------------------
# --- Reaction fluxes and presence ---
# ------------------------------------

# Read FBA data:
fba_normalBlood_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/Normal_Blood.csv',
                          row.names=1, check.names = FALSE)
# Sigmoid Transform fluxes:
fba_normalBlood = 2 * ((1 / (1 + exp(-fba_normalBlood_raw))) - 0.5)

# Reaction presence:
rxn_presence = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/reaction_presence.csv',
                        row.names=1, header=TRUE, check.names=FALSE)



# --------------------------------------------------------------
# --- Reaction fluxes and presence, only reactions with GPRs ---
# --------------------------------------------------------------

# 1. GPRs:
gprs = read.csv('./GENERAL/utility_data/HumanGEM-1.8.0_consistent_GPRs.txt', sep='\t', header=FALSE,
                row.names = 1)
# 1.1. Get reactions with GPRs
rxns_with_gpr = rownames(gprs)[gprs$V2 != ''] # 7 348
# 1.2. From the reactions with GPRs, get only those with at least one gene in the rna dataset
rxns_keep = unique(genes_pathways$reaction[genes_pathways$gene_symbol%in%rownames(cpm_counts_tcells)])
# 1.2.2. 6 881 reactions will be kept

# 2. Get fluxes dataset only with these reactions
fba_normalBlood_gprs = fba_normalBlood[rxns_keep, ]
write.csv(fba_normalBlood_gprs,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/fba_normalBlood_gprs.csv')

# 5. Get reaction presence only with the reactions
rxn_presence_gprs = rxn_presence[rxns_keep, ]
write.csv(rxn_presence_gprs,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/rxn_presence_gprs.csv')



# ------------
# --- UMAPs---
# -------------

# 0. Load data:
tcells_allGenes = SeuratDisk::LoadH5Seurat('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/tcells_allGenes.h5Seurat')
tcells_metabolicGenes = SeuratDisk::LoadH5Seurat('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/tcells_metabolicGenes.h5Seurat')

CPMs = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/CPMs.csv',
                row.names = 1, check.names = FALSE)
GASs = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/GASs.csv',
                row.names = 1, check.names = FALSE)
RASs = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/RASs.csv',
                row.names = 1, check.names = FALSE)

reaction_presence = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/reaction_presence.csv',
                             row.names = 1, check.names = FALSE)
reaction_presence_gprs = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/rxn_presence_gprs.csv',
                                  row.names = 1, check.names = FALSE)

fba_normalBlood = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/Normal_Blood.csv',
                           row.names=1, check.names = FALSE)
fba_normalBlood = 2 * ((1 / (1 + exp(-fba_normalBlood))) - 0.5)
fba_normalBlood_gprs = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/fba_normalBlood_gprs.csv',
                                row.names = 1, check.names = FALSE)

fba_tumourBlood = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/Tumour_Blood.csv',
                           row.names=1, check.names = FALSE)
fba_tumourBlood = 2 * ((1 / (1 + exp(-fba_tumourBlood))) - 0.5)
fba_tumourBlood_gprs = fba_tumourBlood[row.names(fba_normalBlood_gprs), ]


# 1. UMAP of scRNAseq dataset with all genes
Seurat::DimPlot(tcells_allGenes, group.by='models_annotation', label=FALSE, pt.size=.8) +
  ggplot2::scale_colour_manual(values=ct_colors) + ggplot2::ggtitle('scRNAseq | All Genes') +
  ggplot2::theme_light() +
  ggplot2::theme(legend.position = 'bottom', legend.text=ggplot2::element_text(size=8)) +
  ggplot2::guides(col = ggplot2::guide_legend(nrow = 3, override.aes = list(size=6)))


# 2. UMAP of scRNAseq dataset with only the metabolic genes
tcells_metabolicGenes = Seurat::NormalizeData(tcells_metabolicGenes, assay='RNA')
tcells_metabolicGenes = Seurat::ScaleData(tcells_metabolicGenes)
tcells_metabolicGenes = Seurat::RunPCA(tcells_metabolicGenes, assay='integrated')
tcells_metabolicGenes = Seurat::RunUMAP(tcells_metabolicGenes, dims = 1:30, verbose = FALSE)
invisible(gc())
Seurat::DimPlot(tcells_metabolicGenes, group.by='models_annotation', label=FALSE, pt.size=.8) +
  ggplot2::scale_colour_manual(values=ct_colors) + ggplot2::ggtitle('scRNAseq | Metabolic Genes') +
  ggplot2::theme_light() +
  ggplot2::theme(legend.position = 'bottom', legend.text=ggplot2::element_text(size=8)) +
  ggplot2::guides(col = ggplot2::guide_legend(nrow = 3, override.aes = list(size=6)))


# 3. UMAP of CPMs dataset
CPMs_seurat = SeuratObject::CreateSeuratObject(CPMs, meta.data = metadata[colnames(CPMs),])
CPMs_seurat = Seurat::ScaleData(CPMs_seurat)
CPMs_seurat = Seurat::RunPCA(CPMs_seurat, assay='RNA', features=rownames(CPMs_seurat))
CPMs_seurat = Seurat::RunUMAP(CPMs_seurat, dims = 1:30, verbose = FALSE)
invisible(gc())
Seurat::DimPlot(CPMs_seurat, group.by='cell_type', label=FALSE, pt.size=1.8) +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::ggtitle('pseudo-bulk (CPMs)') +
  ggplot2::theme_light() +
  ggplot2::theme(legend.position = 'bottom', legend.text=ggplot2::element_text(size=8)) +
  ggplot2::guides(col = ggplot2::guide_legend(nrow = 3, override.aes = list(size=6)))





# -----------------------------------------------------
# --- Boxplot of number of reactions per cell-types ---
# -----------------------------------------------------

m2 = metadata
m2$cell_type = factor(m2$cell_type, levels = names(ct_colors))
ggplot2::ggplot(m2, ggplot2::aes(x=state, y=n_reactions, colour=cell_type)) +
  ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::facet_wrap(ggplot2::vars(cell_type)) +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::ylab('Number Reactions') + ggplot2::xlab('') +
  ggplot2::theme_light() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1)) #+
  #ggplot2::theme(legend.position = 'none')





# ------------------------
# --- Machine Learning ---
# ------------------------


# 0. Train and test samples:
train_test_samples = list()
train_test_samples$train = sample(colnames(CPMs), 139)
train_test_samples$test = colnames(CPMs)[!colnames(CPMs)%in%train_test_samples$train]
# Distribution of cell-types in train samples
table(metadata[train_test_samples$train, 'cell_type']) / table(metadata[, 'cell_type']) * 100
jsonlite::write_json(train_test_samples,
                     './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/train_test_samples.json')
# Where data will be stored:
ml_res = data.frame(MCC=rep(NA, 7), row.names=c('CPMs', 'GASs', 'RASs',
                                                 'reaction_presence', 'reaction_presence_gprs',
                                                 'fba_normalBlood', 'fba_normalBlood_gprs'))
ml_varImp = list()
# Read it:
train_test_samples = jsonlite::read_json('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/train_test_samples.json',
                                         simplifyVector=TRUE)



# 1. CPMs
# 1.1. Feature Selection
# 1.1.1. Remove low variance reactions
genes_to_remove = caret::nearZeroVar(t(CPMs)[train_test_samples$train,])
# 1.2. Random Forest
# 1.2.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_cpms = caret::train(x = t(CPMs)[train_test_samples$train,-genes_to_remove],
                                y = as.factor(metadata[train_test_samples$train, 'cell_type']),
                                method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_cpms, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/1_cpms_model.RData')
# 1.2.2. Predict test data:
predictions_rf_cpms = predict(cv_model_rf_cpms, t(CPMs)[train_test_samples$test,])
# 1.2.3. Measure prediction capacity
rf_cpms_cm = caret::confusionMatrix(predictions_rf_cpms, as.factor(metadata[train_test_samples$test, 'cell_type']))$table
write.csv(rf_cpms_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/1_cpms_confusionMatrix.csv')
rf_cpms_mcc = yardstick::mcc(data.frame(pred=predictions_rf_cpms,
                                        obs=as.factor(metadata[train_test_samples$test, 'cell_type'])),
                             'obs', 'pred')
ml_res['CPMs', 'MCC'] = rf_cpms_mcc$.estimate
# 1.2.4. Variable Importance:
ml_varImp[['CPMs']] = caret::varImp(cv_model_rf_cpms)$importance


# 2. GASs
# 2.1. Feature Selection
# 2.1.1. Remove low variance reactions
genes_to_remove = caret::nearZeroVar(t(GASs)[train_test_samples$train,])
# 2.2. Random Forest
# 2.2.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_gas = caret::train(x = t(GASs)[train_test_samples$train,-genes_to_remove],
                                y = as.factor(metadata[train_test_samples$train, 'cell_type']),
                                method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_gas, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/2_gas_model.RData')
# 2.2.2. Predict test data:
predictions_rf_gas = predict(cv_model_rf_gas, t(GASs)[train_test_samples$test,])
# 2.2.3. Measure prediction capacity
rf_gas_cm = caret::confusionMatrix(predictions_rf_gas, as.factor(metadata[train_test_samples$test, 'cell_type']))$table
write.csv(rf_gas_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/2_gas_confusionMatrix.csv')
rf_gas_mcc = yardstick::mcc(data.frame(pred=predictions_rf_gas,
                                        obs=as.factor(metadata[train_test_samples$test, 'cell_type'])),
                             'obs', 'pred')
ml_res['GASs', 'MCC'] = rf_gas_mcc$.estimate
# 2.2.4. Variable Importance:
ml_varImp[['GASs']] = caret::varImp(cv_model_rf_gas)$importance


# 3. RASs
# 3.1. Feature Selection
# 3.1.1. Remove low variance reactions
genes_to_remove = caret::nearZeroVar(t(RASs)[train_test_samples$train,])
# 3.2. Random Forest
# 3.2.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_ras = caret::train(x = t(RASs)[train_test_samples$train,-genes_to_remove],
                               y = as.factor(metadata[train_test_samples$train, 'cell_type']),
                               method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_ras, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/3_ras_model.RData')
# 3.2.2. Predict test data:
predictions_rf_ras = predict(cv_model_rf_ras, t(RASs)[train_test_samples$test,])
# 3.2.3. Measure prediction capacity
rf_ras_cm = caret::confusionMatrix(predictions_rf_ras, as.factor(metadata[train_test_samples$test, 'cell_type']))$table
write.csv(rf_ras_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/3_ras_confusionMatrix.csv')
rf_ras_mcc = yardstick::mcc(data.frame(pred=predictions_rf_ras,
                                        obs=as.factor(metadata[train_test_samples$test, 'cell_type'])),
                             'obs', 'pred')
ml_res['RASs', 'MCC'] = rf_ras_mcc$.estimate
# 3.2.4. Variable Importance:
ml_varImp[['RASs']] = caret::varImp(cv_model_rf_ras)$importance


# 4. Reaction Presence
# 4.1. Feature Selection
# 4.1.1. Remove low variance reactions
genes_to_remove = caret::nearZeroVar(t(reaction_presence)[train_test_samples$train,])
# 4.2. Random Forest
# 4.2.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_reaction_presence = caret::train(x = t(reaction_presence)[train_test_samples$train,-genes_to_remove],
                               y = as.factor(metadata[train_test_samples$train, 'cell_type']),
                               method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_reaction_presence, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/4_reaction_presence_model.RData')
# 4.2.2. Predict test data:
predictions_rf_reaction_presence = predict(cv_model_rf_reaction_presence, t(reaction_presence)[train_test_samples$test,])
# 4.2.3. Measure prediction capacity
rf_reaction_presence_cm = caret::confusionMatrix(predictions_rf_reaction_presence, as.factor(metadata[train_test_samples$test, 'cell_type']))$table
write.csv(rf_reaction_presence_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/4_reaction_presence_confusionMatrix.csv')
rf_reaction_presence_mcc = yardstick::mcc(data.frame(pred=predictions_rf_reaction_presence,
                                        obs=as.factor(metadata[train_test_samples$test, 'cell_type'])),
                             'obs', 'pred')
ml_res['reaction_presence', 'MCC'] = rf_reaction_presence_mcc$.estimate
# 4.2.4. Variable Importance:
ml_varImp[['reaction_presence']] = caret::varImp(cv_model_rf_reaction_presence)$importance


# 5. Reaction Presence, only those with GPRs
# 5.1. Feature Selection
# 5.1.1. Remove low variance reactions
genes_to_remove = caret::nearZeroVar(t(reaction_presence_gprs)[train_test_samples$train,])
# 5.2. Random Forest
# 5.2.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_reaction_presence_gprs = caret::train(x = t(reaction_presence_gprs)[train_test_samples$train,-genes_to_remove],
                                             y = as.factor(metadata[train_test_samples$train, 'cell_type']),
                                             method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_reaction_presence_gprs, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/5_reaction_presence_gprs_model.RData')
# 5.2.2. Predict test data:
predictions_rf_reaction_presence_cprs = predict(cv_model_rf_reaction_presence_gprs, t(reaction_presence_gprs)[train_test_samples$test,])
# 5.2.3. Measure prediction capacity
rf_reaction_presence_gprs_cm = caret::confusionMatrix(predictions_rf_reaction_presence_cprs, as.factor(metadata[train_test_samples$test, 'cell_type']))$table
write.csv(rf_reaction_presence_gprs_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/5_reaction_presence_gprs_confusionMatrix.csv')
rf_reaction_presence_gprs_mcc = yardstick::mcc(data.frame(pred=predictions_rf_reaction_presence_cprs,
                                        obs=as.factor(metadata[train_test_samples$test, 'cell_type'])),
                             'obs', 'pred')
ml_res['reaction_presence_gprs', 'MCC'] = rf_reaction_presence_gprs_mcc$.estimate
# 5.2.4. Variable Importance:
ml_varImp[['reaction_presence_gprs']] = caret::varImp(cv_model_rf_reaction_presence_gprs)$importance


# 6. FBA fluxes in normal blood medium
# 6.1. Feature Selection
# 6.1.1. Remove low variance reactions
genes_to_remove = caret::nearZeroVar(t(fba_normalBlood)[train_test_samples$train,])
# 6.2. Random Forest
# 6.2.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_fba_normalBlood = caret::train(x = t(fba_normalBlood)[train_test_samples$train,-genes_to_remove],
                                                  y = as.factor(metadata[train_test_samples$train, 'cell_type']),
                                                  method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_fba_normalBlood, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/6_fba_normalBlood_model.RData')
# 6.2.2. Predict test data:
predictions_rf_fba_normalBlood = predict(cv_model_rf_fba_normalBlood, t(fba_normalBlood)[train_test_samples$test,])
# 6.2.3. Measure prediction capacity
rf_fba_normalBlood_cm = caret::confusionMatrix(predictions_rf_fba_normalBlood, as.factor(metadata[train_test_samples$test, 'cell_type']))$table
write.csv(rf_fba_normalBlood_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/6_fba_normalBlood_confusionMatrix.csv')
rf_fba_normalBlood_mcc = yardstick::mcc(data.frame(pred=predictions_rf_fba_normalBlood,
                                        obs=as.factor(metadata[train_test_samples$test, 'cell_type'])),
                             'obs', 'pred')
ml_res['fba_normalBlood', 'MCC'] = rf_fba_normalBlood_mcc$.estimate
# 6.2.4. Variable Importance:
ml_varImp[['fba_normalBlood']] = caret::varImp(cv_model_rf_fba_normalBlood)$importance


# 7. FBA fluxes in normal blood medium, only those with GPRs
# 7.1. Feature Selection
# 7.1.1. Remove low variance reactions
genes_to_remove = caret::nearZeroVar(t(fba_normalBlood_gprs)[train_test_samples$train,])
# 7.2. Random Forest
# 7.2.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_fba_normalBlood_gprs = caret::train(x = t(fba_normalBlood_gprs)[train_test_samples$train,-genes_to_remove],
                                           y = as.factor(metadata[train_test_samples$train, 'cell_type']),
                                           method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_fba_normalBlood_gprs, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/7_fba_normalBlood_gprs_model.RData')
# 7.2.2. Predict test data:
predictions_rf_fba_normalBlood_gprs = predict(cv_model_rf_fba_normalBlood_gprs, t(fba_normalBlood_gprs)[train_test_samples$test,])
# 7.2.3. Measure prediction capacity
rf_fba_normalBlood_gprs_cm = caret::confusionMatrix(predictions_rf_fba_normalBlood_gprs, as.factor(metadata[train_test_samples$test, 'cell_type']))$table
write.csv(rf_fba_normalBlood_gprs_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/7_fba_normalBlood_gprs_confusionMatrix.csv')
rf_fba_normalBlood_gprs_mcc = yardstick::mcc(data.frame(pred=predictions_rf_fba_normalBlood_gprs,
                                        obs=as.factor(metadata[train_test_samples$test, 'cell_type'])),
                             'obs', 'pred')
ml_res['fba_normalBlood_gprs', 'MCC'] = rf_fba_normalBlood_gprs_mcc$.estimate
# 7.2.4. Variable Importance:
ml_varImp[['fba_normalBlood_gprs']] = caret::varImp(cv_model_rf_fba_normalBlood_gprs)$importance


# 8. FBA fluxes in normal blood medium
# 8.1. Feature Selection
# 8.1.1. Remove low variance reactions
genes_to_remove = caret::nearZeroVar(t(fba_tumourBlood)[train_test_samples$train,])
# 8.2. Random Forest
# 8.2.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_fba_tumourBlood = caret::train(x = t(fba_normalBlood)[train_test_samples$train,-genes_to_remove],
                                           y = as.factor(metadata[train_test_samples$train, 'cell_type']),
                                           method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_fba_tumourBlood, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/6_fba_tumourBlood_model.RData')
# 8.2.2. Predict test data:
predictions_rf_fba_tumourBlood = predict(cv_model_rf_fba_tumourBlood, t(fba_normalBlood)[train_test_samples$test,])
# 8.2.3. Measure prediction capacity
rf_fba_tumourBlood_cm = caret::confusionMatrix(predictions_rf_fba_tumourBlood, as.factor(metadata[train_test_samples$test, 'cell_type']))$table
write.csv(rf_fba_tumourBlood_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/6_fba_tumourBlood_confusionMatrix.csv')
rf_fba_tumourBlood_mcc = yardstick::mcc(data.frame(pred=predictions_rf_fba_tumourBlood,
                                                   obs=as.factor(metadata[train_test_samples$test, 'cell_type'])),
                                        'obs', 'pred')
ml_res['fba_tumourBlood', 'MCC'] = rf_fba_tumourBlood_mcc$.estimate
# 8.2.4. Variable Importance:
ml_varImp[['fba_tumourBlood']] = caret::varImp(cv_model_rf_fba_tumourBlood)$importance


# 9. FBA fluxes in normal blood medium, only those with GPRs
# 9.1. Feature Selection
# 9.1.1. Remove low variance reactions
genes_to_remove = caret::nearZeroVar(t(fba_tumourBlood_gprs)[train_test_samples$train,])
# 9.2. Random Forest
# 9.2.1. Train:
fitControl = caret::trainControl(method="repeatedcv", number=10, repeats=10)
cv_model_rf_fba_tumourBlood_gprs = caret::train(x = t(fba_tumourBlood_gprs)[train_test_samples$train,-genes_to_remove],
                                                y = as.factor(metadata[train_test_samples$train, 'cell_type']),
                                                method = 'rf', trControl = fitControl)
saveRDS(cv_model_rf_fba_tumourBlood_gprs, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/7_fba_tumourBlood_gprs_model.RData')
# 9.2.2. Predict test data:
predictions_rf_fba_tumourBlood_gprs = predict(cv_model_rf_fba_tumourBlood_gprs, t(fba_tumourBlood_gprs)[train_test_samples$test,])
# 9.2.3. Measure prediction capacity
rf_fba_tumourBlood_gprs_cm = caret::confusionMatrix(predictions_rf_fba_tumourBlood_gprs, as.factor(metadata[train_test_samples$test, 'cell_type']))$table
write.csv(rf_fba_tumourBlood_gprs_cm, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/7_fba_tumourBlood_gprs_confusionMatrix.csv')
rf_fba_tumourBlood_gprs_mcc = yardstick::mcc(data.frame(pred=predictions_rf_fba_tumourBlood_gprs,
                                                        obs=as.factor(metadata[train_test_samples$test, 'cell_type'])),
                                             'obs', 'pred')
ml_res['fba_tumourBlood_gprs', 'MCC'] = rf_fba_tumourBlood_gprs_mcc$.estimate
# 9.2.4. Variable Importance:
ml_varImp[['fba_tumourBlood_gprs']] = caret::varImp(cv_model_rf_fba_tumourBlood_gprs)$importance


# 10. Store MCC metric values and variable importance list:
write.csv(ml_res, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/MCC.csv')
jsonlite::write_json(ml_varImp, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/ml/varImp.json')


# 11. Visualize MCCs
plot_df = data.frame(mcc=ml_res$MCC, dataset=rownames(ml_res))
# 11.1 Just with CPMs, FBA normal blood and reaction presence (GPRs)
plot_df2 = plot_df[c(1,5,7),]
plot_df2$dataset = factor(c('Pseudo-Bulk RNA (CPMs)', 'Reaction Presence', 'FBA Normal Blood'),
                          levels=c('Pseudo-Bulk RNA (CPMs)', 'Reaction Presence', 'FBA Normal Blood'))
ggplot2::ggplot(plot_df2, ggplot2::aes(y=mcc, x=dataset)) +
  ggplot2::geom_bar(stat="identity", fill='#E69F00', colour='black', width=0.5) +
  ggplot2::theme_minimal() + ggplot2::xlab('') + ggplot2::ylab('MCC') + ggplot2::ylim(-.05, 1) #+
  #ggplot2::theme(axis.text.x=ggplot2::element_text(angle=20, hjust=1, vjust=1))
