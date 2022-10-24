ct_colors = c('#b30000', '#999999', '#cc79a7', '#b3c9e6', '#4e84c4', '#c3d7a4', '#52854c', '#fdb981',
              '#d16103', '#293352')
names(ct_colors) = c('Cytotoxic CD8 Tcells', 'Follicular CD4 Tcells', 'IL17+ CD4 Tcells',
                     'Memory CD4 Tcells', 'Memory CD8 Tcells', 'Naive CD4 Tcells', 'Naive CD8 Tcells',
                     'Proliferative CD4 Tcells', 'Proliferative CD8 Tcells',  'Regulatory CD4 Tcells')
state_colors = c('#52854c', '#b30000')
names(state_colors) = c('Normal Matched', 'Tumour')
cms_colors = c('#d4cfcf', '#cc79a7', '#b3c9e6', '#c3d7a4', '#52854c', '#fdb981')
names(cms_colors) = c('CMS3', 'Mixed', 'CMS1', 'CMS2', 'Normal Matched', 'CMS4')



# --------------------------------------
# --- Load Gene Essentiality Results ---
# --------------------------------------

#gene_essentiality = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Gene_essentiality/growth.csv',
#                             row.names=1, header=TRUE, check.names=FALSE)
gene_essentiality = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Gene_essentiality/growth_biomass.csv',
                             row.names=1, header=TRUE, check.names=FALSE)
# Change to gene symbols:
humangem_genemapping = unlist(jsonlite::read_json('./GENERAL/utility_data/genes_mapping.json',
                                                  simplifyVector=TRUE))
rownames(gene_essentiality) = names(humangem_genemapping)[match(rownames(gene_essentiality),
                                                                humangem_genemapping)]

# Also, load original FBA, under normal conditions:
fba_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/Normal_Blood_biomass.csv',
                          row.names=1, header=TRUE, check.names = FALSE)[,colnames(gene_essentiality)]

# Metadata:
metadata = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/metadata.csv', row.names=1)

# House-keeping genes:
hk_genes = read.csv('./GENERAL/utility_data/Housekeeping_GenesHuman.csv')



# ------------------------------------------------------
# --- Get list of essential genes for each cell-type ---
# ------------------------------------------------------

# 1. Genes that are essential for each model (objective turns completely zero):
# 1.1. Create matrix where results will be stored:
model_essentiality = matrix(rep(FALSE, dim(gene_essentiality)[1]*dim(gene_essentiality)[2]),
                            nrow=dim(gene_essentiality)[1])
rownames(model_essentiality) = rownames(gene_essentiality)
colnames(model_essentiality) = colnames(gene_essentiality)
# 1.2. Essential genes in a model will be TRUE, others will be FALSE:
for(gene in rownames(model_essentiality)){
  message(gene)
  outcome = gene_essentiality[gene, colnames(model_essentiality)]
  outcome2 = (outcome <= 0 | is.na(outcome))
  model_essentiality[gene, ] = outcome2
}
# 1.3. Remove genes known to be housekeeping genes:
model_essentiality_woHK = model_essentiality[!rownames(model_essentiality)%in%hk_genes$Gene.name, ]
#write.csv(model_essentiality_woHK,
#          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Gene_essentiality/model_essential.csv')
write.csv(model_essentiality_woHK,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Gene_essentiality/model_essential_biomass.csv')
# 1.4. Genes essential for each cell-type:
cellType_essentiality_woHK = list()
for(ct in unique(metadata[, 'cell_type'])){
  message(ct)
  models_names = rownames(metadata)[metadata$cell_type==ct]
  n_models = length(models_names)
  ct_essentiality = rowSums(model_essentiality_woHK[, models_names])
  cellType_essentiality_woHK[[ct]] = rownames(model_essentiality_woHK)[(ct_essentiality / n_models) > .5]
}
#jsonlite::write_json(cellType_essentiality_woHK,
#                     './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Gene_essentiality/CT_essential.json')
jsonlite::write_json(cellType_essentiality_woHK,
                     './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Gene_essentiality/CT_essential_biomass.json')





# -----------------------------
# --- Map genes to pathways ---
# -----------------------------

# 1. Get gene -> reaction dataframe:
# 1.1. Read GPRs file:
gprs = read.csv('./GENERAL/utility_data/HumanGEM-1.8.0_consistent_GPRs.txt', sep='\t',
                header=FALSE, row.names=1)
# 1.2. Remove reactions with no GPR:
gprs_genes = rownames(gprs)[gprs$V2!='']
gprs = gprs[gprs$V2!='','V2']
names(gprs) = gprs_genes
# 1.3. Create dataframe where columns will be gene and reaction:
genes = c()
rxns = c()
for(rxn in names(gprs)){
  gpr = gprs[rxn]
  gpr = gsub('[()]', '', gpr)
  gpr_split = unlist(strsplit(strsplit(gpr, ' or ')[[1]], ' and '))
  gpr_split = gsub(' ', '', gpr_split)
  genes = c(genes, gpr_split)
  rxns = c(rxns, rep(rxn, length(gpr_split)))
}
genes_symbols = names(humangem_genemapping)[match(genes, humangem_genemapping)]
genes_pathways = data.frame(gene=genes, gene_symbol=genes_symbols, reaction=rxns)
# 1.4. Add mapping between reaction and pathway:
rxns_subsystems = jsonlite::read_json('./GENERAL/utility_data/reactions_subsystems_mapping.json',
                                      simplifyVector=TRUE)
genes_pathways$pathways = NA
for(i in 1:dim(genes_pathways)[1]){
  rxn = genes_pathways$reaction[i]
  if(!is.null(rxns_subsystems[[rxn]]))
    genes_pathways$pathways[i] = rxns_subsystems[[rxn]]
}
genes_pathways = genes_pathways[!is.na(genes_pathways$pathways),]
# 1.5. Save file:
write.csv(genes_pathways, './GENERAL/utility_data/genes_subsystems_mapping.csv', row.names=FALSE)



# -------------------------------------------------------------
# --- Percentage of genes from a pathway that are essential ---
# -------------------------------------------------------------

# 0.
genes_pathways = read.csv('./GENERAL/utility_data/genes_subsystems_mapping.csv')

# 1. Get only reactions that are affected when the gene is deleted:
single_gene_deletions = jsonlite::read_json('./GENERAL/utility_data/single_gene_deletion.json',
                                            simplifyVector = TRUE)
reactions_affected = unique(unlist(single_gene_deletions))
genes_pathways_affected = genes_pathways[genes_pathways$reaction%in%reactions_affected,]

# 2. For the population of those reactions, see the most affected pathways. From all the genes in
#    a pathway whose deletion leads to inactivation of at least one reaction in the pathway, what is
#    the percentage of essential genes
n_genes_pathway = table(genes_pathways_affected$pathways)
res_mat = matrix(rep(0, length(unique(genes_pathways_affected$pathways)) * dim(model_essentiality_woHK)[2]),
                 ncol=dim(model_essentiality_woHK)[2])
colnames(res_mat) = colnames(model_essentiality_woHK)
rownames(res_mat) = unique(genes_pathways_affected$pathways)
res_mat_n = res_mat
for(path in unique(genes_pathways_affected$pathways)){
  message(path)
  for(model in colnames(model_essentiality_woHK)){
    n_gns = sum(rownames(model_essentiality_woHK)[model_essentiality_woHK[, model]] %in% 
                  genes_pathways_affected$gene_symbol[genes_pathways_affected$pathways==path])
    perc = n_gns / n_genes_pathway[path] * 100
    res_mat[path, model] = perc
    res_mat_n[path, model] = n_gns
  }
}
# Save file:
#write.csv(res_mat,
#          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Gene_essentiality/percentage_essential_pathways_models.csv')
write.csv(res_mat,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Gene_essentiality/percentage_essential_pathways_models_biomass.csv')
write.csv(res_mat_n,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Gene_essentiality/ngenes_essential_pathways_models_biomass.csv')

# 3. Heatmap of most affected pathways:
ht3 = ComplexHeatmap::Heatmap(res_mat[, rownames(metadata)],
                              col = circlize::colorRamp2(c(0, 50, 100),
                                                         c('blue', 'white', 'red')),
                              name='% Essential Genes',
                              column_split = metadata$cell_type, column_title = NULL,
                              show_column_names = FALSE,
                              row_names_gp = grid::gpar(fontsize=7.5),
                              row_names_max_width = ComplexHeatmap::max_text_width(rownames(res_mat), 
                                                                                   gp = grid::gpar(fontsize = 12)),
                              heatmap_legend_param = list(direction='horizontal'),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata$cell_type, state=metadata$state, cms=metadata$CMS, col=list(cell_type=ct_colors, state=state_colors, cms=cms_colors), annotation_legend_param = list(cell_type=list(title='Cell Type', ncol=5), state=list(title='State', nrow=1), cms=list(title='CMS', nrow=1)), show_annotation_name=FALSE, simple_anno_size = grid::unit(2,'mm')))
ComplexHeatmap::draw(ht3, merge_legend = TRUE, heatmap_legend_side='bottom',
                     annotation_legend_side='bottom')

medians = matrixStats::rowMedians(as.matrix(res_mat))
names(medians) = rownames(res_mat)
top_20 = sort(medians, decreasing=T)[1:20]
ht3 = ComplexHeatmap::Heatmap(res_mat[names(top_20), rownames(metadata)],
                              col = circlize::colorRamp2(c(0, 50, 100),
                                                         c('blue', 'white', 'red')),
                              name='% Essential Genes',
                              column_split = metadata$cell_type, column_title = NULL,
                              show_column_names = FALSE,
                              row_names_gp = grid::gpar(fontsize=7.5),
                              row_names_max_width = ComplexHeatmap::max_text_width(rownames(res_mat), 
                                                                                   gp = grid::gpar(fontsize = 12)),
                              heatmap_legend_param = list(direction='horizontal'),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata$cell_type, state=metadata$state, cms=metadata$CMS, col=list(cell_type=ct_colors, state=state_colors, cms=cms_colors), annotation_legend_param = list(cell_type=list(title='Cell Type', ncol=5), state=list(title='State', nrow=1), cms=list(title='CMS', nrow=1)), show_annotation_name=FALSE, simple_anno_size = grid::unit(2,'mm')))
ComplexHeatmap::draw(ht3, merge_legend = TRUE, heatmap_legend_side='bottom',
                     annotation_legend_side='bottom')

# Nº essential genes per pathway
ht3 = ComplexHeatmap::Heatmap(res_mat_n[, rownames(metadata)],
                              #col = circlize::colorRamp2(c(0, 50, 100),
                              #                           c('blue', 'white', 'red')),
                              name='N Essential Genes',
                              column_split = metadata$cell_type, column_title = NULL,
                              show_column_names = FALSE,
                              row_names_gp = grid::gpar(fontsize=7.5),
                              row_names_max_width = ComplexHeatmap::max_text_width(rownames(res_mat), 
                                                                                   gp = grid::gpar(fontsize = 12)),
                              heatmap_legend_param = list(direction='horizontal'),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata$cell_type, state=metadata$state, cms=metadata$CMS, col=list(cell_type=ct_colors, state=state_colors, cms=cms_colors), annotation_legend_param = list(cell_type=list(title='Cell Type', ncol=5), state=list(title='State', nrow=1), cms=list(title='CMS', nrow=1)), show_annotation_name=FALSE, simple_anno_size = grid::unit(2,'mm')))
ComplexHeatmap::draw(ht3, merge_legend = TRUE, heatmap_legend_side='bottom',
                     annotation_legend_side='bottom')

medians = matrixStats::rowMedians(as.matrix(res_mat_n))
names(medians) = rownames(res_mat_n)
top_20_n = sort(medians, decreasing=T)[1:21]
ht3 = ComplexHeatmap::Heatmap(res_mat_n[names(top_20_n), rownames(metadata)],
                              #col = circlize::colorRamp2(c(0, 50, 100),
                              #                           c('blue', 'white', 'red')),
                              name='N Essential Genes',
                              column_split = metadata$cell_type, column_title = NULL,
                              show_column_names = FALSE,
                              row_names_gp = grid::gpar(fontsize=7.5),
                              row_names_max_width = ComplexHeatmap::max_text_width(rownames(res_mat), 
                                                                                   gp = grid::gpar(fontsize = 12)),
                              heatmap_legend_param = list(direction='horizontal'),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata$cell_type, state=metadata$state, cms=metadata$CMS, col=list(cell_type=ct_colors, state=state_colors, cms=cms_colors), annotation_legend_param = list(cell_type=list(title='Cell Type', ncol=5), state=list(title='State', nrow=1), cms=list(title='CMS', nrow=1)), show_annotation_name=FALSE, simple_anno_size = grid::unit(2,'mm')))
ComplexHeatmap::draw(ht3, merge_legend = TRUE, heatmap_legend_side='bottom',
                     annotation_legend_side='bottom')
ht3 = ComplexHeatmap::Heatmap(res_mat_n[names(top_20_n[2:21]), rownames(metadata)],
                              #col = circlize::colorRamp2(c(0, 50, 100),
                              #                           c('blue', 'white', 'red')),
                              name='N Essential Genes',
                              column_split = metadata$cell_type, column_title = NULL,
                              show_column_names = FALSE,
                              row_names_gp = grid::gpar(fontsize=7.5),
                              row_names_max_width = ComplexHeatmap::max_text_width(rownames(res_mat), 
                                                                                   gp = grid::gpar(fontsize = 12)),
                              heatmap_legend_param = list(direction='horizontal'),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata$cell_type, state=metadata$state, cms=metadata$CMS, col=list(cell_type=ct_colors, state=state_colors, cms=cms_colors), annotation_legend_param = list(cell_type=list(title='Cell Type', ncol=5), state=list(title='State', nrow=1), cms=list(title='CMS', nrow=1)), show_annotation_name=FALSE, simple_anno_size = grid::unit(2,'mm')))
ComplexHeatmap::draw(ht3, merge_legend = TRUE, heatmap_legend_side='bottom',
                     annotation_legend_side='bottom')

# Essential genes that are transport reactions:
transport_genes = unique(genes_pathways$gene_symbol[genes_pathways$pathways=='Transport reactions'])
essential_genes_transports = cellType_essentiality_woHK
for(ct in names(essential_genes_transports)){
  essential_genes_transports[[ct]] = essential_genes_transports[[ct]][essential_genes_transports[[ct]]%in%transport_genes]
}
sort(unique(unlist(essential_genes_transports)))





# ----------------------------------------------------------------------
# --- Representation of essential genes from a pathway in the models ---
# ----------------------------------------------------------------------

# 1. Get data:
res_mat2 = matrix(rep(0, length(unique(genes_pathways_affected$pathways)) * dim(model_essentiality_woHK)[2]),
                 ncol=dim(model_essentiality_woHK)[2])
colnames(res_mat2) = colnames(model_essentiality_woHK)
rownames(res_mat2) = unique(genes_pathways_affected$pathways)
for(path in unique(genes_pathways_affected$pathways)){
  message(path)
  for(model in colnames(model_essentiality_woHK)){
    n_total_essential_genes = sum(model_essentiality_woHK[,model])
    n_gns = sum(rownames(model_essentiality_woHK)[model_essentiality_woHK[, model]] %in% 
                  genes_pathways_affected$gene_symbol[genes_pathways_affected$pathways==path])
    perc = n_gns / n_total_essential_genes * 100
    res_mat2[path, model] = perc
  }
}
# Save file:
#write.csv(res_mat2,
#          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Gene_essentiality/representation_essential_pathways_models.csv')
write.csv(res_mat2,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Gene_essentiality/representation_essential_pathways_models_biomass.csv')

# 2. Heatmap:
ht3 = ComplexHeatmap::Heatmap(res_mat2[, rownames(metadata)],
                              col = circlize::colorRamp2(c(0, max(res_mat2)/2,
                                                           max(res_mat2)),
                                                         c('blue', 'white', 'red')),
                              name='% Presence',
                              column_split = metadata$cell_type,
                              column_title = NULL,
                              show_column_names = FALSE,
                              row_names_gp = grid::gpar(fontsize=7.5),
                              heatmap_legend_param = list(direction='horizontal'),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata$cell_type, state=metadata$state, cms=metadata$CMS, col=list(cell_type=ct_colors, state=state_colors, cms=cms_colors), annotation_legend_param = list(cell_type=list(title='Cell Type', ncol=5), state=list(title='State', nrow=1), cms=list(title='CMS', nrow=1)), show_annotation_name=FALSE, simple_anno_size = grid::unit(2,'mm')))
ComplexHeatmap::draw(ht3, merge_legend = TRUE, heatmap_legend_side='bottom',
                     annotation_legend_side='bottom')


# Heatmap of potencial essential genes from certain pathways:
genes = c('FECH', 'PPOX', 'COQ3', 'COQ6', 'COQ7', 'RFK', 'BLVRB',
          'ALOX12B', 'CYP2F1', 'CYP4F8', 'CYP4F12', 'LTC4S', 'HPGD',
          "A4GALT", "GLA", "B3GALNT1", "B3GALT5", "GBGT1", "NAGA", "ST8SIA1", "ACAT2", "AACS", "HMGCS1", "ECHDC2",
          "ALB", "SERPINA3", "SERPINA1", "APOB", "APOC1", "APOC2", "APOC3", "FGA", "PLG", "F2", "TFRC", "APOE", "BMP1",
          "CARNS1", "HNMT", "HAL", "FTCD", "AMDHD1", "UROC1",
          'SLC12A3', 'SLC6A20', 'SLC6A4', 'SLC29A4', 'SLC7A10', 'SLC7A5', 'SLC38A5', 'SLC1A4', 'SLC7A11',
          'SLC29A2', 'SLC28A3', 'SLC28A1', 'SLC5A6')
genes_paths = c(rep('Heme synthesis', 2), rep('Ubiquinone synthesis', 3), rep('Riboflavin metabolism', 2),
                rep('Eicosanoid metabolism', 6),
                rep('Glycosphingolipid biosynthesis-globo series', 7), rep('Butanoate metabolism', 4),
                rep('Protein assembly', 13),
                rep('Histidine metabolism', 6),
                rep('Chloride uptake', 2), rep('Serotonin uptake', 2), rep('Uptake of amino acids', 5),
                rep('Nucleotides/nucleosides uptake', 3), 'Vitamins B5 and B7, and lipoic acid uptake')
genes_map = model_essentiality_woHK[genes, rownames(metadata)]
genes_map[genes_map==TRUE] = 1
hte = ComplexHeatmap::Heatmap(genes_map,
                              show_heatmap_legend = FALSE,
                              col = circlize::colorRamp2(c(0, 1),
                                                         c('blue', 'red')),
                              name='Essential',
                              column_split = metadata$cell_type, row_split = genes_paths,
                              column_title = NULL, row_title = NULL,
                              show_column_names = FALSE,
                              cluster_rows = TRUE,
                              #row_names_gp = grid::gpar(fontsize=5.5),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata$cell_type, state=metadata$state, cms=metadata$CMS, col=list(cell_type=ct_colors, state=state_colors, cms=cms_colors), annotation_legend_param = list(cell_type=list(title='Cell Type', ncol=5), state=list(title='State', nrow=1), cms=list(title='CMS', nrow=1)), show_annotation_name=FALSE, simple_anno_size = grid::unit(2,'mm')))
ComplexHeatmap::draw(hte, merge_legend = TRUE, heatmap_legend_side='bottom',
                     annotation_legend_side='bottom')



# ----------------------
# --- Uniquely genes ---
# ----------------------

# 1. Get genes uniquely essential in each cell-type:
all_cts = names(cellType_essentiality_woHK)
unique_essential = list()
for(ct in all_cts){
  ct_essential_genes = cellType_essentiality_woHK[[ct]]
  other_ct_essential_genes = unique(unlist(cellType_essentiality_woHK[all_cts[all_cts!=ct]]))
  ct_unique_essential_genes = ct_essential_genes[!ct_essential_genes%in%other_ct_essential_genes]
  unique_essential[[ct]] = ct_unique_essential_genes
}
#jsonlite::write_json(unique_essential,
#                     './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Gene_essentiality/unique_essential_ct.json')
jsonlite::write_json(unique_essential,
                     './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Gene_essentiality/unique_essential_ct_biomass.json')

# 2. Get genes uniquely decrease in each cell-type:
all_cts = names(cellType_decrease_woHK)
unique_decrease = list()
for(ct in all_cts){
  ct_essential_genes = cellType_decrease_woHK[[ct]]
  other_ct_essential_genes = unique(unlist(cellType_decrease_woHK[all_cts[all_cts!=ct]]))
  ct_unique_essential_genes = ct_essential_genes[!ct_essential_genes%in%other_ct_essential_genes]
  unique_decrease[[ct]] = ct_unique_essential_genes
}
# No uniques





# -------------------------------------------------
# --- Comparison with gene essentiality studies ---
# -------------------------------------------------

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
cd8_cells = c('Cytotoxic CD8 Tcells', 'Proliferative CD8 Tcells',
              'Memory CD8 Tcells', 'Naive CD8 Tcells')
cd4_cells = c('Proliferative CD4 Tcells', 'Memory CD4 Tcells', 'Naive CD4 Tcells',
              'Regulatory CD4 Tcells', 'Follicular CD4 Tcells', 'IL17+ CD4 Tcells')



# ---
# - Read results obtained above
# ---

# House-keeping genes
hk_genes = read.csv('./GENERAL/utility_data/Housekeeping_GenesHuman.csv')
# Essential genes by cell-type
#cellType_essentiality_woHK = jsonlite::read_json('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Gene_essentiality/CT_essential.json', simplifyVector=TRUE)
cellType_essentiality_woHK = jsonlite::read_json('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Gene_essentiality/CT_essential_biomass.json', simplifyVector=TRUE)
# Essential genes - CD8
cd8_insilico = unique(unlist(cellType_essentiality_woHK[cd8_cells]))
# Essential genes - CD4
cd4_insilico = unique(unlist(cellType_essentiality_woHK[cd4_cells]))
# Gene-pathways mapping
genes_pathways = read.csv('./GENERAL/utility_data/genes_subsystems_mapping.csv')
# Genes tested
single_genes_dels = jsonlite::read_json('./GENERAL/utility_data/single_gene_deletion.json', simplifyVector = T)
gp2 = unique(genes_pathways[,c('gene', 'gene_symbol')])
genes_tested = gp2$gene_symbol[gp2$gene%in%names(single_genes_dels)]
# The genes tested that are house keeping genes
hk_genes_tested = genes_tested[genes_tested%in%hk_genes$Gene.name]
# Genes tested without HK genes
genes_tested_woHK = genes_tested[!genes_tested%in%hk_genes_tested]



# ---
# - Read literature Data
# ---

# Shifrut et al - human CD8
essential_genes_cd8_shifrut_r = 
  xlsx::read.xlsx('./GENERAL/utility_data/gene_essential_studies/CD8_human_shifrut_etal.xlsx',
                  1)
essential_genes_cd8_shifrut_r = essential_genes_cd8_shifrut_r[!is.na(essential_genes_cd8_shifrut_r$id),]

all_genes_cd8lit = essential_genes_cd8_shifrut_r$id

essential_genes_cd8_shifrut_r = essential_genes_cd8_shifrut_r[essential_genes_cd8_shifrut_r$neg.lfc <= (-0.1),]
essential_genes_cd8_shifrut_genes = unique(essential_genes_cd8_shifrut_r$id)

essential_genes_cd8_shifrut_metabolic = essential_genes_cd8_shifrut_genes[essential_genes_cd8_shifrut_genes
                                                                          %in%genes_pathways$gene_symbol]
essential_genes_cd8_shifrut_metabolic_hk = essential_genes_cd8_shifrut_metabolic[essential_genes_cd8_shifrut_metabolic%in%hk_genes$Gene.name]
essential_genes_cd8_shifrut_metabolic_woHK = essential_genes_cd8_shifrut_metabolic[!essential_genes_cd8_shifrut_metabolic%in%essential_genes_cd8_shifrut_metabolic_hk]


# Ting et al - human CD4
essential_genes_cd4_ting_r = 
  xlsx::read.xlsx('./GENERAL/utility_data/gene_essential_studies/cd4_CRISPR_tingetal.xlsx', 1,
                  startRow=3)
all_genes_cd4lit = essential_genes_cd4_ting_r$Gene_ID

essential_genes_cd4_ting_r = essential_genes_cd4_ting_r[!is.na(essential_genes_cd4_ting_r$X.log2FC..max.),]
essential_genes_cd4_ting_r = essential_genes_cd4_ting_r[essential_genes_cd4_ting_r$X.logRSA>3 &
                                                          essential_genes_cd4_ting_r$X.log2FC..max.>0,]

essential_genes_cd4_ting_genes = unique(essential_genes_cd4_ting_r$Gene_ID)
essential_genes_cd4_ting_metabolic = essential_genes_cd4_ting_genes[essential_genes_cd4_ting_genes
                                                                    %in%genes_pathways$gene_symbol]

essential_genes_cd4_ting_metabolic_hk = essential_genes_cd4_ting_metabolic[essential_genes_cd4_ting_metabolic%in%hk_genes$Gene.name]
essential_genes_cd4_ting_metabolic_woHK = essential_genes_cd4_ting_metabolic[!essential_genes_cd4_ting_metabolic%in%essential_genes_cd4_ting_metabolic_hk]



# ---
# - Genes tested vs tested in lit CD4 vs tested in lit CD8
# ---

display_venn(list(cd4litall=all_genes_cd4lit, cd8litall=all_genes_cd8lit,
                  genes_tested=genes_tested),
             category.names = c("CD4 (Ting et al, 2018)" , "CD8 (Shifrut et al, 2018)", "Genes Tested"),
             # Circles
             fill = c("darkblue", "darkgreen", "grey"),
             # Numbers
             cex = 1.5,
             fontface = "italic",
             # Set names
             cat.cex = 1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(-0.04, -0.03, -0.37),
             cat.pos=c(-10, 7, 0),
             # Title
             main='',
             main.pos = c(.5, .85),
             main.cex = 1.5,
             main.fontface = 'bold')



# ---
# - Genes tested vs literature CD4 vs literature CD8
# ---

display_venn(list(cd4lit=essential_genes_cd4_ting_metabolic_woHK,
                  cd8slit=essential_genes_cd8_shifrut_metabolic_woHK,
                  genes_tested=genes_tested),
             category.names = c("CD4 (Ting et al, 2018)" , "CD8 (Shifrut et al, 2018)", "Genes Tested"),
             # Circles
             fill = c("#56B4E9", "#009E73", "grey"),
             # Numbers
             cex = 1.5,
             fontface = "italic",
             # Set names
             cat.cex = 1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(-0.04, -0.04, -.37),
             cat.pos=c(-10, 10, 0),
             # Title
             main='',
             main.pos = c(.5, .85),
             main.cex = 1.5,
             main.fontface = 'bold')



# ---
# - Genes tested vs silico CD4 vs silico CD8
# ---

display_venn(list(cd4silico=cd4_insilico, cd8silico=cd8_insilico,
                  genes_tested=genes_tested),
             category.names = c("CD4 Predictions" , "CD8 Predictions", "Genes Tested"),
             # Circles
             fill = c("#E69F00", "#990033", "grey"),
             # Numbers
             cex = 1.5,
             fontface = "italic",
             # Set names
             cat.cex = 1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(-0.04, -0.04, -0.04),
             cat.pos=c(-10, 0, 10),
             # Title
             main='',
             main.pos = c(.5, .85),
             main.cex = 1.5,
             main.fontface = 'bold')



# ---
# - silico CD4 vs literature CD4 (common genes tested)
# ---

cd4lit_tested = essential_genes_cd4_ting_metabolic_woHK[essential_genes_cd4_ting_metabolic_woHK%in%
                                                          genes_tested_woHK]
cd4pred_inlit = cd4_insilico[cd4_insilico%in%all_genes_cd4lit]
display_venn(list(cd4lit=cd4lit_tested,
                  cd4silico=cd4pred_inlit),
             category.names = c("CD4 (Ting et al, 2018)" , "CD4 Predictions"),
             ext.percent=.01,
             # Circles
             fill = c("#56B4E9", "#E69F00"),
             # Numbers
             cex = 1.5,
             fontface = "italic",
             # Set names
             cat.cex = 1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(-.21, -.415),
             cat.pos=c(-10, 0),
             # Title
             main='',
             main.pos = c(.5, .85),
             main.cex = 1.5,
             main.fontface = 'bold')



# ---
# - silico CD8 vs literature CD8 (common genes tested)
# ---

cd8lit_tested = essential_genes_cd8_shifrut_metabolic_woHK[essential_genes_cd8_shifrut_metabolic_woHK%in%
                                                             genes_tested_woHK]
cd8pred_inlit = cd8_insilico[cd8_insilico%in%all_genes_cd8lit]
display_venn(list(cd8lit=cd8lit_tested,
                  cd8silico=cd8pred_inlit),
             category.names = c("CD8 (Shifrut et al, 2018)" , "CD8 Predictions"),
             ext.percent=.1,
             # Circles
             fill = c("#009E73", "#990033"),
             # Numbers
             cex = 1.5,
             fontface = "italic",
             # Set names
             cat.cex = 1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(.012, .012),
             cat.pos=c(0, 0),
             # Title
             main='',
             main.pos = c(.5, .85),
             main.cex = 1.5,
             main.fontface = 'bold')
#



# ---
# - literature CD4 and CD8 genes (from common genes tested) not predicted by us: do they still decrease?
# ---

cd4_notpredicted = cd4lit_tested[!cd4lit_tested%in%cd4pred_inlit]
cd8_notpredicted = cd8lit_tested[!cd8lit_tested%in%cd8pred_inlit]

cd8_cells = c('Cytotoxic CD8 Tcells', 'Proliferative CD8 Tcells',
              'Memory CD8 Tcells', 'Naive CD8 Tcells')
cd4_cells = c('Proliferative CD4 Tcells', 'Memory CD4 Tcells', 'Naive CD4 Tcells',
              'Regulatory CD4 Tcells', 'Follicular CD4 Tcells', 'IL17+ CD4 Tcells')

# CD4:
cd4_zeros = c()
cd4_decreases = c()
cd4_cts = c()
for(gene in cd4_notpredicted){
  for(ct in cd4_cells){
    models_names = grep(ct, colnames(fba_fluxes_raw), value=T)
    cd4_cts = c(cd4_cts, ct)
    
    z = as.numeric(gene_essentiality[gene, models_names]) == 0
    z[is.na(z)] = TRUE
    cd4_zeros = c(cd4_zeros, sum(z) / 196 * 100)
    
    d = as.numeric(gene_essentiality[gene, models_names]) < as.numeric(fba_fluxes_raw['MAR13082', models_names])
    d[is.na(d)] = TRUE
    cd4_decreases = c(cd4_decreases, sum(d) / 196 * 100)
  }
  
}
df_plot = data.frame(gene=cd4_notpredicted, cell_type=cd4_cts, cd4_zeros, cd4_decreases)

# CD8:
cd8_zeros = c()
cd8_decreases = c()
cd8_cts = c()
for(gene in cd8_notpredicted){
  for(ct in cd8_cells){
    models_names = grep(ct, colnames(fba_fluxes_raw), value=T)
    cd8_cts = c(cd8_cts, ct)
    
    z = as.numeric(gene_essentiality[gene, models_names]) == 0
    z[is.na(z)] = TRUE
    cd8_zeros = c(cd8_zeros, sum(z) / 196 * 100)
    
    d = as.numeric(gene_essentiality[gene, models_names]) < as.numeric(fba_fluxes_raw['MAR13082', models_names])
    d[is.na(d)] = TRUE
    cd8_decreases = c(cd8_decreases, sum(d) / 196 * 100)
  }
  
}
df_plot_cd8 = data.frame(gene=cd8_notpredicted, cell_type=cd8_cts, cd8_zeros, cd8_decreases)





# ----------------------------------------------
# --- Essential genes of transport reactions ---
# ----------------------------------------------

# Potential essential:
transport = unique(genes_pathways_affected$gene_symbol[genes_pathways_affected$pathways=='Transport reactions'])

# Essential transports: 43
transport_essential = list()
for(ct in names(cellType_essentiality_woHK)){
  genes = cellType_essentiality_woHK[[ct]]
  transport_essential[[ct]] = genes[genes%in%transport]
}
x = unique(unlist(transport_essential))

g = 'SLC5A8'
g%in%essential_genes_cd4_ting_metabolic_woHK
g%in%all_genes_cd4lit

g%in%essential_genes_cd8_shifrut_metabolic_woHK
g%in%all_genes_cd8lit

unlist(transport_essential)[unlist(transport_essential)=='SLC7A11']



