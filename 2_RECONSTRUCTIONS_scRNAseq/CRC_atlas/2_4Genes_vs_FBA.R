ct_colors = c('#b30000', '#999999', '#cc79a7', '#b3c9e6', '#4e84c4', '#c3d7a4', '#52854c', '#fdb981',
              '#d16103', '#293352')
names(ct_colors) = c('Cytotoxic CD8 Tcells', 'Follicular CD4 Tcells', 'IL17+ CD4 Tcells',
                     'Memory CD4 Tcells', 'Memory CD8 Tcells', 'Naive CD4 Tcells', 'Naive CD8 Tcells',
                     'Proliferative CD4 Tcells', 'Proliferative CD8 Tcells',  'Regulatory CD4 Tcells')

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

pseudo_bulk_tcells = NULL
tcelltypes_metabolic = readRDS('./0Data/scRNAseq/CRC_atlas/expression_data.Rdata')
for(indiv in names(tcelltypes_metabolic)){
  if(indiv%in%names(models_created$NormalMatched)){
    for(samp in names(tcelltypes_metabolic[[indiv]])){
      message(samp)
      expression = tcelltypes_metabolic[[indiv]][[samp]]
      expression = expression[, colnames(expression)%in%models_created$NormalMatched[[indiv]]$control[[samp]]]
      expression = cbind(rownames(expression), expression)
      colnames(expression) = c('ID', paste(indiv, samp, colnames(expression)[-1], sep='_'))
      if(is.null(pseudo_bulk_tcells)) pseudo_bulk_tcells = expression
      else pseudo_bulk_tcells = merge(pseudo_bulk_tcells, expression, by='ID')
    }
  }
}
rownames(pseudo_bulk_tcells) = pseudo_bulk_tcells$ID
pseudo_bulk_tcells = pseudo_bulk_tcells[metabolic_genes[metabolic_genes%in%rownames(pseudo_bulk_tcells)],-1]
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

# 3. Get reaction activity only with these reactions
fba_activity_gprs = fba_fluxes_abs_gprs
fba_activity_gprs[fba_activity_gprs>0] = 1
fba_activity_gprs[fba_activity_gprs<0] = -1

# 4. Get reaction presence only with the reactions
rxn_presence = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/general/reaction_presence_Tcells.csv',
                        row.names=1, header=TRUE, check.names=FALSE)
rxn_presence_gprs = rxn_presence[rxns_keep, ]





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



