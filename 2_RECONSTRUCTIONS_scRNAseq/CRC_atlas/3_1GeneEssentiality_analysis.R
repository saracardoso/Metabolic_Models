ct_colors = c('#b30000', '#999999', '#cc79a7', '#b3c9e6', '#4e84c4', '#c3d7a4', '#52854c', '#fdb981',
              '#d16103', '#293352')
names(ct_colors) = c('Cytotoxic CD8 Tcells', 'Follicular CD4 Tcells', 'IL17+ CD4 Tcells',
                     'Memory CD4 Tcells', 'Memory CD8 Tcells', 'Naive CD4 Tcells', 'Naive CD8 Tcells',
                     'Proliferative CD4 Tcells', 'Proliferative CD8 Tcells',  'Regulatory CD4 Tcells')



# --------------------------------------
# --- Load Gene Essentiality Results ---
# --------------------------------------

gene_essentiality = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/3_gene_essentiality/origObj_growth.csv',
                             row.names=1, header=TRUE, check.names=FALSE)
# Change to gene symbols:
humangem_genemapping = unlist(jsonlite::read_json('./GENERAL/utility_data/genes_mapping.json',
                                                  simplifyVector=TRUE))
rownames(gene_essentiality) = names(humangem_genemapping)[match(rownames(gene_essentiality),
                                                                humangem_genemapping)]

# Also, load original FBA, under normal conditions:
fba_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/normal_FBA.csv',
                          row.names=1, header=TRUE, check.names = FALSE)[,colnames(gene_essentiality)]

# Metadata:
split_colnames = strsplit(colnames(fba_fluxes_raw), '_')
samples = c()
cell_types = c()
individuals = c()
for(row_split in split_colnames){
  individuals = c(individuals, row_split[1])
  samples = c(samples, row_split[2])
  cell_types = c(cell_types, row_split[3])
}
metadata = data.frame(individual=individuals, sample=samples, cell_type=cell_types,
                      row.names=colnames(fba_fluxes_raw))
more_meta = read.csv('./0Data/scRNAseq/CRC_atlas/expression_data/metadata.csv', row.names=1)
states = c()
for(samp in metadata$sample){
  states = c(states, more_meta[samp, 'Sample.Source'])
}
metadata$state = states

# House-keeping genes:
hk_genes = read.csv('./GENERAL/utility_data/Housekeeping_GenesHuman.csv')



# ------------------------------------------------------
# --- Get list of essential genes for each cell-type ---
# ------------------------------------------------------

# 0. No models need to be removed from analysis, as no proliferative models had original biomass zero
# and no other models had biomass + atp production zero.
prol_tcells = grep('Proliferative', colnames(fba_fluxes_raw), value=TRUE)
other_tcells = grep('Proliferative', colnames(fba_fluxes_raw), value=TRUE, invert=TRUE)
sum(fba_fluxes_raw['MAR13082', prol_tcells] <= 0)
sum(colSums(fba_fluxes_raw[c('MAR13082', 'MAR06916'), prol_tcells]) <= 0)

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
write.csv(model_essentiality_woHK,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/3_gene_essentiality/model_essential.csv')
# 1.4. Genes essential for each cell-type:
cellType_essentiality_woHK = list()
for(ct in unique(metadata[, 'cell_type'])){
  message(ct)
  models_names = rownames(metadata)[metadata$cell_type==ct]
  n_models = length(models_names)
  ct_essentiality = rowSums(model_essentiality_woHK[, models_names])
  cellType_essentiality_woHK[[ct]] = rownames(model_essentiality_woHK)[ct_essentiality / n_models > .9]
}
jsonlite::write_json(cellType_essentiality_woHK,
                     './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/3_gene_essentiality/CT_essential.json')

# 2. Genes that are not essential for each model (objective turns completely zero) but still decrease enough:
# 2.1. Create matrix where results will be stored:
model_decrease= matrix(rep(FALSE, dim(gene_essentiality)[1]*dim(gene_essentiality)[2]),
                            nrow=dim(gene_essentiality)[1])
rownames(model_decrease) = rownames(gene_essentiality)
colnames(model_decrease) = colnames(gene_essentiality)
# 2.2. Essential genes in a model will be TRUE, others will be FALSE:
for(gene in rownames(model_decrease)){
  message(gene)
  outcome = gene_essentiality[gene, colnames(model_decrease)]
  original = rep(0, length(colnames(model_decrease)))
  names(original) = colnames(model_decrease)
  original[prol_tcells] = fba_fluxes_raw['MAR13082', prol_tcells]
  original[other_tcells] = colSums(fba_fluxes_raw[c('MAR13082', 'MAR06916'), other_tcells])
  outcome2 = (outcome < as.numeric(original)/2 | is.na(outcome))
  model_decrease[gene, ] = outcome2
}
# 2.3. Remove genes known to be housekeeping genes:
model_decrease_woHK = model_decrease[!rownames(model_decrease)%in%hk_genes$Gene.name, ]
write.csv(model_decrease_woHK,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/3_gene_essentiality/model_decrease.csv')
# 1.4. Genes essential for each cell-type (only test genes not in ct's essential):
cellType_decrease_woHK = list()
for(ct in unique(metadata[, 'cell_type'])){
  message(ct)
  models_names = rownames(metadata)[metadata$cell_type==ct]
  n_models = length(models_names)
  ct_essentiality = rowSums(model_decrease_woHK[, models_names])
  cellType_decrease_woHK[[ct]] = rownames(model_decrease_woHK)[ct_essentiality / n_models > .9]
  cellType_decrease_woHK[[ct]] = cellType_decrease_woHK[[ct]][!cellType_decrease_woHK[[ct]]%in%
                                                                cellType_essentiality_woHK[[ct]]]
}
jsonlite::write_json(cellType_decrease_woHK,
                     './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/3_gene_essentiality/CT_decrease.json')


# 3. Genes whose objective increases:
# 3.1. Create matrix where results will be stored:
model_increase = matrix(rep(FALSE, dim(gene_essentiality)[1]*dim(gene_essentiality)[2]),
                        nrow=dim(gene_essentiality)[1])
rownames(model_increase) = rownames(gene_essentiality)
colnames(model_increase) = colnames(gene_essentiality)
# 3.2. Increasing genes in models:
for(gene in rownames(model_increase)){
  message(gene)
  outcome = gene_essentiality[gene, ]
  original = rep(0, length(colnames(model_increase)))
  names(original) = colnames(model_increase)
  original[prol_tcells] = fba_fluxes_raw['MAR13082', prol_tcells]
  original[other_tcells] = colSums(fba_fluxes_raw[c('MAR13082', 'MAR06916'), other_tcells])
  outcome2 = (outcome > 0 & !is.na(outcome) & outcome > original)
  model_increase[gene, ] = outcome2
}
write.csv(model_increase,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/3_gene_essentiality/model_increase.csv')
# 3.3. Increasing genes in cell-types:
cellType_increase = list()
for(ct in unique(metadata$cell_type)){
  message(ct)
  models_names = rownames(metadata)[metadata$cell_type==ct]
  n_models = length(models_names)
  ct_essentiality = rowSums(model_increase[, models_names])
  cellType_increase[[ct]] = rownames(model_increase)[ct_essentiality / n_models > 0.9]
}
# No increases.



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
for(path in unique(genes_pathways_affected$pathways)){
  message(path)
  for(model in colnames(model_essentiality_woHK)){
    n_gns = sum(rownames(model_essentiality_woHK)[model_essentiality_woHK[, model]] %in% 
                  genes_pathways_affected$gene_symbol[genes_pathways_affected$pathways==path])
    perc = n_gns / n_genes_pathway[path] * 100
    res_mat[path, model] = perc
  }
}
# Save file:
write.csv(res_mat,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/3_gene_essentiality/percentage_essential_pathways_models.csv')

# 3. Heatmap of most affected pathways:
stt = metadata$state
stt[stt!='Normal Matched'] = 'Tumour'
state_colors2 = c('#52854c', '#b30000')
names(state_colors2) = c('Normal Matched', 'Tumour')
ht3 = ComplexHeatmap::Heatmap(t(res_mat),
                              col = circlize::colorRamp2(c(0, 50, 100),
                                                         c('blue', 'white', 'red')),
                              name='% Essential Genes',
                              row_split = metadata$cell_type,
                              row_title = NULL,
                              show_row_names = FALSE,
                              column_names_gp = grid::gpar(fontsize=7.5),
                              column_names_rot = 65,
                              heatmap_legend_param = list(direction='horizontal'),
                              left_annotation = ComplexHeatmap::rowAnnotation(cell_type=metadata$cell_type, state=stt, col=list(cell_type=ct_colors, state=state_colors2), annotation_legend_param = list(cell_type=list(title='Cell Type', ncol=5), state=list(title='State', nrow=1)), show_annotation_name=FALSE, simple_anno_size = grid::unit(2,'mm')))

ComplexHeatmap::draw(ht3, merge_legend = TRUE, heatmap_legend_side='bottom',
                     annotation_legend_side='bottom')



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
write.csv(res_mat2,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/3_gene_essentiality/representation_essential_pathways_models.csv')

# 2. Heatmap:
stt = metadata$state
stt[stt!='Normal Matched'] = 'Tumour'
state_colors2 = c('#52854c', '#b30000')
names(state_colors2) = c('Normal Matched', 'Tumour')
ht3 = ComplexHeatmap::Heatmap(res_mat2,
                              col = circlize::colorRamp2(c(0, max(res_mat2)/2,
                                                           max(res_mat2)),
                                                         c('blue', 'white', 'red')),
                              name='% Presence',
                              column_split = metadata$cell_type,
                              column_title = NULL,
                              show_column_names = FALSE,
                              row_names_gp = grid::gpar(fontsize=5.5),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata$cell_type, state=stt, col=list(cell_type=ct_colors, state=state_colors2), annotation_legend_param = list(cell_type=list(title='Cell Type'), state=list(title='State')), show_annotation_name=FALSE, simple_anno_size = grid::unit(2,'mm')))
ComplexHeatmap::draw(ht3, merge_legend = TRUE, heatmap_legend_side='bottom',
                     annotation_legend_side='bottom')

# Heatmap of potencial essential genes from eicosanoid metabolism
genes = genes_pathways_affected$gene_symbol[genes_pathways_affected$pathways=='Eicosanoid metabolism']
stt = metadata$state
stt[stt!='Normal Matched'] = 'Tumour'
state_colors2 = c('#52854c', '#b30000')
names(state_colors2) = c('Normal Matched', 'Tumour')
genes_map = model_essentiality_woHK[genes, rownames(metadata)]
genes_map[genes_map==TRUE] = 1
hte = ComplexHeatmap::Heatmap(genes_map,
                              show_heatmap_legend = FALSE,
                              col = circlize::colorRamp2(c(0, 1),
                                                         c('lightgrey', 'red')),
                              name='Essential',
                              column_split = metadata$cell_type,
                              column_title = NULL,
                              show_column_names = FALSE,
                              cluster_rows = TRUE,
                              #row_names_gp = grid::gpar(fontsize=5.5),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata$cell_type,
                                                                                 state=stt, col=list(cell_type=ct_colors, state=state_colors2),
                                                                                 annotation_legend_param = list(cell_type=list(title='Cell Type'),
                                                                                                                state=list(title='State')),
                                                                                 show_annotation_name=FALSE, simple_anno_size = grid::unit(2,'mm')))
ComplexHeatmap::draw(hte, merge_legend = TRUE, #heatmap_legend_side='none',
                     annotation_legend_side='bottom')

# Heatmap of potencial essential genes from eicosanoid metabolism, leukotriene metabolism, prostaglandin biosynthesis
genes = unique(genes_pathways_affected$gene_symbol[genes_pathways_affected$pathways%in%
                                                     c('Eicosanoid metabolism',
                                                       'Leukotriene metabolism',
                                                       'Prostaglandin biosynthesis')])
stt = metadata$state
stt[stt!='Normal Matched'] = 'Tumour'
state_colors2 = c('#52854c', '#b30000')
names(state_colors2) = c('Normal Matched', 'Tumour')
genes_map = model_essentiality_woHK[genes[-3], rownames(metadata)]
genes_map[genes_map==TRUE] = 1
hte = ComplexHeatmap::Heatmap(genes_map,
                              show_heatmap_legend = FALSE,
                              col = circlize::colorRamp2(c(0, 1),
                                                         c('lightgrey', 'red')),
                              name='Essential',
                              column_split = metadata$cell_type,
                              column_title = NULL,
                              show_column_names = FALSE,
                              cluster_rows = TRUE,
                              #row_names_gp = grid::gpar(fontsize=5.5),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata$cell_type,
                                                                                 state=stt, col=list(cell_type=ct_colors, state=state_colors2),
                                                                                 annotation_legend_param = list(cell_type=list(title='Cell Type'),
                                                                                                                state=list(title='State')),
                                                                                 show_annotation_name=FALSE, simple_anno_size = grid::unit(2,'mm')))
ComplexHeatmap::draw(hte, merge_legend = TRUE, #heatmap_legend_side='none',
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
jsonlite::write_json(unique_essential,
                     './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/3_gene_essentiality/unique_essential_ct.json')

# 2. Get genes uniquely decrease in each cell-type:
all_cts = names(cellType_decrease_woHK)
unique_decrease = list()
for(ct in all_cts){
  ct_essential_genes = cellType_decrease_woHK[[ct]]
  other_ct_essential_genes = unique(unlist(cellType_decrease_woHK[all_cts[all_cts!=ct]]))
  ct_unique_essential_genes = ct_essential_genes[!ct_essential_genes%in%other_ct_essential_genes]
  unique_decrease[[ct]] = ct_unique_essential_genes
}
jsonlite::write_json(unique_decrease,
                     './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/3_gene_essentiality/unique_decrease_ct.json')





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
cellType_essentiality_woHK = jsonlite::read_json('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/3_gene_essentiality/CT_essential.json', simplifyVector=TRUE)
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
# - Genes tested vs silico CD4 vs silico CD8
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
# - silico CD4 vs literature CD4 (- not tested)
# ---

cd4lit_tested = essential_genes_cd4_ting_metabolic_woHK[essential_genes_cd4_ting_metabolic_woHK%in%
                                                          genes_tested_woHK]
display_venn(list(cd4lit=cd4lit_tested,
                  cd4silico=cd4_insilico),
             category.names = c("CD4 (Ting et al, 2018)" , "CD4 Predictions"),
             ext.percent=.001,
             # Circles
             fill = c("#56B4E9", "#E69F00"),
             # Numbers
             cex = 1.5,
             fontface = "italic",
             # Set names
             cat.cex = 1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(0.02, -0.37),
             cat.pos=c(40, 0),
             # Title
             main='',
             main.pos = c(.5, .85),
             main.cex = 1.5,
             main.fontface = 'bold')

# Were the genes predicted as essential but not essential in study studied (or has info for it)?
pred_cd4_notinLit = cd4_insilico[!cd4_insilico%in%cd4lit_tested]
sum(pred_cd4_notinLit%in%all_genes_cd4lit) # 174

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
             cat.dist = c(-.125, -.415),
             cat.pos=c(-10, 0),
             # Title
             main='',
             main.pos = c(.5, .85),
             main.cex = 1.5,
             main.fontface = 'bold')



# ---
# - silico CD8 vs literature CD8 (- not tested)
# ---

cd8lit_tested = essential_genes_cd8_shifrut_metabolic_woHK[essential_genes_cd8_shifrut_metabolic_woHK%in%
                                                             genes_tested_woHK]
display_venn(list(cd8lit=cd8lit_tested,
                  cd8silico=cd8_insilico),
             category.names = c("CD8 (Shifrut et al, 2018)" , "CD8 Predictions"),
             ext.percent=.001,
             # Circles
             fill = c("#009E73", "#990033"),
             # Numbers
             cex = 1.5,
             fontface = "italic",
             # Set names
             cat.cex = 1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(0.01, -0.37),
             cat.pos=c(0, 0),
             # Title
             main='',
             main.pos = c(.5, .85),
             main.cex = 1.5,
             main.fontface = 'bold')

# Were the genes predicted as essential but not essential in study studied (or has info for it)?
pred_cd8_notinLit = cd8_insilico[!cd8_insilico%in%cd8lit_tested]
sum(pred_cd8_notinLit%in%all_genes_cd8lit) # 126

cd8pred_inlit = cd8_insilico[cd8_insilico%in%all_genes_cd8lit]
display_venn(list(cd8lit=cd8lit_tested,
                  cd8silico=cd8pred_inlit),
             category.names = c("Essential" , "CD8 Predictions"),
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
             cat.dist = c(-.37, -.42),
             cat.pos=c(0, 0),
             # Title
             main='',
             main.pos = c(.5, .85),
             main.cex = 1.5,
             main.fontface = 'bold')
#





# -------------------------------------------------------------
# --- Uniquely essential genes - representation in pathways ---
# -------------------------------------------------------------



# ---
# - Read results obtained above
# ---

# Essential genes by cell-type
cellType_essentiality_woHK = jsonlite::read_json('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/3_gene_essentiality/CT_essential.json', simplifyVector=TRUE)
# Gene-pathways mapping
genes_pathways = read.csv('./GENERAL/utility_data/genes_subsystems_mapping.csv')
# Genes uniquely essential for each cell-type
unique_essential = jsonlite::read_json('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/3_gene_essentiality/unique_essential_ct.json', simplifyVector = T)


# ---
# - Get number of genes and reactions in pathways affected by unique essential genes for each ct
# ---

cts = c()
paths = c()
ngenes = c()
nrxns = c()
for(ct in names(unique_essential)){
  x = genes_pathways[genes_pathways$gene_symbol%in%cellType_essentiality_woHK[[ct]],]
  x$unique = NA
  x[x$gene_symbol%in%unique_essential[[ct]],'unique'] = 'yes'
  x = x[!is.na(x$unique), ]
  
  pathways = unique(x$pathways)
  for(path in pathways){
    y = x[x$pathways==path, ]
    n_genes = length(unique(y$gene_symbol))
    n_rxns = length(unique(y$reaction))
    
    cts = c(cts, ct)
    paths = c(paths, path)
    ngenes = c(ngenes, n_genes)
    nrxns = c(nrxns, n_rxns)
  }
}
df_plot = data.frame(cell_type=cts, pathway=paths, n_gene=ngenes, n_reaction=nrxns, ratio=nrxns/ngenes)


# ---
# - Plot the information
# ---

tab_paths = table(df_plot$pathway)

ggplot2::ggplot(df_plot[df_plot$pathway%in%names(tab_paths[tab_paths<3]),],
                ggplot2::aes(pathway, cell_type)) +
  ggplot2::geom_point(ggplot2::aes(size=n_gene, color=ratio)) +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1)) +
  ggplot2::scale_color_continuous(type='viridis')


