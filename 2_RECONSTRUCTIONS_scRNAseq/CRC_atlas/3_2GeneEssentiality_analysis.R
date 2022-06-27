ct_colors = c('#b30000', '#999999', '#cc79a7', '#b3c9e6', '#4e84c4', '#c3d7a4', '#52854c', '#fdb981',
              '#d16103', '#293352')
names(ct_colors) = c('Cytotoxic CD8 Tcells', 'Follicular CD4 Tcells', 'IL17+ CD4 Tcells',
                     'Memory CD4 Tcells', 'Memory CD8 Tcells', 'Naive CD4 Tcells', 'Naive CD8 Tcells',
                     'Proliferative CD4 Tcells', 'Proliferative CD8 Tcells',  'Regulatory CD4 Tcells')
state_colors = c('#52854c', '#b30000')
names(state_colors) = c('Normal Matched', 'Tumour')





# --------------------------------------------------------------------------------
# --- COMPARE PSEUDO-BULK, TAS, RAS AND REACTION PRESENCE FOR CERTAIN PATHWAYS ---
# --------------------------------------------------------------------------------


# ---
# - Get the data
# ---

# 0. Metabolic genes:
genes_pathways = read.csv('./GENERAL/utility_data/genes_subsystems_mapping.csv')

# 1. Get the data:
# 1.1. Reaction presence
rxn_presence = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/general/reaction_presence_Tcells.csv',
                        row.names = 1, check.names = FALSE)
# 1.2. Bulk CPM counts, TAS and RAS scores
all_expression_data = readRDS('./0Data/scRNAseq/CRC_atlas/expression_data.Rdata')
data_folder = './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched'
cpm_counts_tcells = c()
tas_counts_tcells = c()
ras_counts_tcells = c()
for(indiv in c('31', '32', '33', '35', 'KUL01', 'KUL19', 'KUL21', 'SMC01', 'SMC04', 'SMC06', 'SMC07',
               'SMC08', 'SMC10')){
  message(indiv)
  for(samp in names(all_expression_data[[indiv]])){
    cpm_counts = read.csv(paste(data_folder, '/', indiv, '/data_info/', samp, '_CPM.csv', sep=''),
                          row.names=1, check.names = FALSE)
    tas_counts = read.csv(paste(data_folder, '/', indiv, '/data_info/', samp, '_TAS.csv', sep=''),
                          row.names=1, check.names = FALSE)
    ras_counts = read.csv(paste(data_folder, '/', indiv, '/data_info/', samp, '_RAS.csv', sep=''),
                          row.names=1, check.names = FALSE)
    
    colnames(cpm_counts) = paste(indiv, samp, colnames(cpm_counts), sep='_')
    colnames(tas_counts) = paste(indiv, samp, colnames(tas_counts), sep='_')
    colnames(ras_counts) = paste(indiv, samp, colnames(ras_counts), sep='_')
    
    cpm_counts = cpm_counts[, colnames(cpm_counts)%in%colnames(rxn_presence)]
    tas_counts = tas_counts[, colnames(tas_counts)%in%colnames(rxn_presence)]
    ras_counts = ras_counts[, colnames(ras_counts)%in%colnames(rxn_presence)]
    
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
}
cpm_counts_tcells[is.na(cpm_counts_tcells)] = 0

genes_map = unique(genes_pathways[,c('gene', 'gene_symbol')])
genes_with_info_ensemble = genes_map[genes_map$gene_symbol%in%rownames(cpm_counts_tcells), ]

cpm_counts_tcells = cpm_counts_tcells[genes_with_info_ensemble$gene_symbol,]
tas_counts_tcells = tas_counts_tcells[genes_with_info_ensemble$gene, ]
rownames(tas_counts_tcells) = genes_with_info_ensemble$gene_symbol

# 1.4. FBA fluxes
fba_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/normal_FBA.csv', row.names=1, check.names = FALSE)
fba_fluxes = 2 * ((1 / (1 + exp(-fba_fluxes_raw))) - 0.5)
fba_fluxes = fba_fluxes[, colnames(cpm_counts_tcells)]


# 1.5. Metadata
split_colnames = strsplit(colnames(cpm_counts_tcells), '_')
samples = c()
cell_types = c()
individuals = c()
for(row_split in split_colnames){
  individuals = c(individuals, row_split[1])
  samples = c(samples, row_split[2])
  cell_types = c(cell_types, row_split[3])
}
metadata = data.frame(individual=individuals, sample=samples, cell_type=cell_types,
                      row.names=colnames(cpm_counts_tcells))
more_meta = read.csv('./0Data/scRNAseq/CRC_atlas/expression_data/metadata.csv', row.names=1)
states = c()
for(samp in metadata$sample){
  states = c(states, more_meta[samp, 'Sample.Source'])
}
states[states!='Normal Matched'] = 'Tumour'
metadata$state = states



# ---
# - Eicosanoid metabolism + Leukotriene metabolism + Prostaglandin biosynthesis
# ---

genes = c('PTGS1', 'PTGS2', 'PTGIS', 'TBXAS1', 'PTGDS', 'HPGDS', 'PTGES', 'PTGES2', 'PTGES3', 'AKR1C3',
          'ALOX5', 'LTA4H', 'RNPEP', 'LTC4S', 'GGT6', 'GGT1', 'GGT2', 'GGT7', 'GGT5', 'DPEP1', 'DPEP2',
          'DPEP3')
rxns = c('MAR01305', 'MAR01307', 'MAR08557', 'MAR08558', 'MAR01308', 'MAR01313', 'MAR01315', 'MAR01330', 'MAR08559',
         'MAR01332', 'MAR01333', 'MAR01434', 'MAR01336', 'MAR01338', 'MAR01312', 'MAR08560', 'MAR01329',
         'MAR00958', 'MAR00959', 'MAR01080', 'MAR08550', 'MAR01084', 'MAR08548', 'MAR01085', 'MAR08552')
pathways_genes = c(rep('Prostaglandins', 10), rep('Leukotrienes', 12))
pathways_rxns = c(rep('Prostaglandins', 17), rep('Leukotrienes', 8))

# 1. CPM Counts
# 1.1. With PTGES3
max_cpm = max(cpm_counts_tcells)
ht_genes = ComplexHeatmap::Heatmap(cpm_counts_tcells[genes, ],
                                   col = circlize::colorRamp2(c(0, max_cpm/2, max_cpm),
                                                              c('lightgrey', 'lightgrey', 'red')),
                        show_column_names = FALSE,
                        column_split = metadata[colnames(cpm_counts_tcells),'cell_type'],
                        column_title = NULL,
                        name = 'CPM counts (left heatmap)',
                        cluster_rows = FALSE,
                        left_annotation = ComplexHeatmap::rowAnnotation(Pathways = pathways_genes,
                                                                        col=list(Pathways=c('Leukotrienes'='#FBC02D', 'Prostaglandins'='#E91E63')),
                                                                        show_annotation_name=FALSE, simple_anno_size=grid::unit(2,'mm')),
                        top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata[colnames(cpm_counts_tcells),'cell_type'],
                                                                           state=metadata[colnames(cpm_counts_tcells),'state'], col=list(cell_type=ct_colors, state=state_colors),
                                                                           annotation_legend_param=list(cell_type=list(title='Cell Type'),
                                                                                                        state=list(title='State')),
                                                                           show_annotation_name=FALSE, simple_anno_size=grid::unit(2,'mm')))

# 2. TAS Scores
min_tas = min(tas_counts_tcells[genes, ])
max_tas = max(tas_counts_tcells[genes, ])
ht_tas = ComplexHeatmap::Heatmap(tas_counts_tcells[genes, ],
                        col = circlize::colorRamp2(c(min_tas, 0, max_tas),
                                                   c('lightgrey', 'lightgrey', 'red')),
                        show_column_names = FALSE,
                        column_split = factor(metadata[colnames(cpm_counts_tcells),'cell_type'],
                                              levels=names(ComplexHeatmap::column_order(ht_genes))),
                        column_title = NULL, 
                        show_heatmap_legend = FALSE,
                        cluster_rows = FALSE,
                        column_order = unlist(ComplexHeatmap::column_order(ht_genes)),
                        top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata[colnames(cpm_counts_tcells),'cell_type'],
                                                                           state=metadata[colnames(cpm_counts_tcells),'state'], col=list(cell_type=ct_colors, state=state_colors),
                                                                           annotation_legend_param=list(cell_type=list(title='Cell Type'),
                                                                                                        state=list(title='State')),
                                                                           show_annotation_name=FALSE, simple_anno_size=grid::unit(2,'mm')))

# 3. See both heatmaps side by side
hts = ht_genes + ht_tas
ComplexHeatmap::draw(hts, merge_legend=TRUE)

# 4. RAS Scores
min_ras = min(ras_counts_tcells[rxns, ])
max_ras = max(ras_counts_tcells[rxns, ])
ht_ras = ComplexHeatmap::Heatmap(ras_counts_tcells[rxns, ],
                                 col = circlize::colorRamp2(c(min_ras, 0, max_ras),
                                                            c('lightgrey', 'lightgrey', 'red')),
                                 show_column_names = FALSE,
                                 column_split = factor(metadata[colnames(ras_counts_tcells),'cell_type'],
                                                       levels=names(ComplexHeatmap::column_order(ht_genes))),
                                 column_title = 'RAS', 
                                 show_heatmap_legend = FALSE,
                                 cluster_rows = FALSE,
                                 column_order = unlist(ComplexHeatmap::column_order(ht_genes)),
                                 top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata[colnames(ras_counts_tcells),'cell_type'],
                                                                                    state=metadata[colnames(ras_counts_tcells),'state'], col=list(cell_type=ct_colors, state=state_colors),
                                                                                    annotation_legend_param=list(cell_type=list(title='Cell Type', ncol=5),
                                                                                                                 state=list(title='State', ncol=2)),
                                                                                    show_annotation_name=FALSE, simple_anno_size=grid::unit(2,'mm')),
                                 left_annotation = ComplexHeatmap::rowAnnotation(Pathways = pathways_rxns,
                                                                                 col=list(Pathways=c('Leukotrienes'='#FBC02D', 'Prostaglandins'='#E91E63')),
                                                                                 annotation_legend_param=list(Pathways=list(ncol=2)),
                                                                                 show_annotation_name=FALSE, simple_anno_size=grid::unit(2,'mm')))

# 5. Reaction Presence
ht_rp = ComplexHeatmap::Heatmap(rxn_presence[rxns, ],
                                 col = circlize::colorRamp2(c(0, 1),
                                                            c('lightgrey', 'red')),
                                 show_column_names = FALSE,
                                 column_split = factor(metadata[colnames(ras_counts_tcells),'cell_type'],
                                                       levels=names(ComplexHeatmap::column_order(ht_genes))),
                                 column_title = 'Presence', 
                                 show_heatmap_legend = FALSE,
                                cluster_rows = FALSE,
                                column_order = unlist(ComplexHeatmap::column_order(ht_genes)),
                                 top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata[colnames(ras_counts_tcells),'cell_type'],
                                                                                    state=metadata[colnames(ras_counts_tcells),'state'], col=list(cell_type=ct_colors, state=state_colors),
                                                                                    annotation_legend_param=list(cell_type=list(title='Cell Type', ncol=5),
                                                                                                                 state=list(title='State', ncol=2)),
                                                                                    show_annotation_name=FALSE, simple_anno_size=grid::unit(2,'mm')))

# 6. Reaction Fluxes
ht_flx = ComplexHeatmap::Heatmap(fba_fluxes[rxns, ],
                                col = circlize::colorRamp2(c(-1, 0, 1),
                                                           c('blue', 'white', 'red')),
                                show_column_names = FALSE,
                                column_split = factor(metadata[colnames(ras_counts_tcells),'cell_type'],
                                                      levels=names(ComplexHeatmap::column_order(ht_genes))),
                                column_title = 'Fluxes', 
                                show_heatmap_legend = FALSE,
                                cluster_rows = FALSE,
                                column_order = unlist(ComplexHeatmap::column_order(ht_genes)),
                                top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata[colnames(ras_counts_tcells),'cell_type'],
                                                                                   state=metadata[colnames(ras_counts_tcells),'state'], col=list(cell_type=ct_colors, state=state_colors),
                                                                                   annotation_legend_param=list(cell_type=list(title='Cell Type', ncol=5),
                                                                                                                state=list(title='State', ncol=2)),
                                                                                   show_annotation_name=FALSE, simple_anno_size=grid::unit(2,'mm')))

#6. Reactions with Fluxes
xx = fba_fluxes[rxns, ]
xx[xx<0] = -1
xx[xx>0] = 1
ht_flx_p = ComplexHeatmap::Heatmap(xx,
                                 col = circlize::colorRamp2(c(-1, 0, 1),
                                                            c('blue', 'lightgrey', 'red')),
                                 show_column_names = FALSE,
                                 column_split = factor(metadata[colnames(ras_counts_tcells),'cell_type'],
                                                       levels=names(ComplexHeatmap::column_order(ht_genes))),
                                 column_title = 'Reactions with Flux', 
                                 show_heatmap_legend = FALSE,
                                 cluster_rows = FALSE,
                                 column_order = unlist(ComplexHeatmap::column_order(ht_genes)),
                                 top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata[colnames(ras_counts_tcells),'cell_type'],
                                                                                    state=metadata[colnames(ras_counts_tcells),'state'], col=list(cell_type=ct_colors, state=state_colors),
                                                                                    annotation_legend_param=list(cell_type=list(title='Cell Type', ncol=5),
                                                                                                                 state=list(title='State', ncol=2)),
                                                                                    show_annotation_name=FALSE, simple_anno_size=grid::unit(2,'mm')))


# 7. See both reaction heatmaps side by side
hts = ht_ras + ht_rp + ht_flx + ht_flx_p
ComplexHeatmap::draw(hts, merge_legend=TRUE, heatmap_legend_side='bottom')


# 8. Boxplots

df = data.frame(values=as.numeric(cpm_counts_tcells['PTGS2', ]), models=colnames(cpm_counts_tcells),
                cell_types=metadata[colnames(cpm_counts_tcells), 'cell_type'],
                state=metadata[colnames(cpm_counts_tcells), 'state'])

ggplot2::ggplot(df, ggplot2::aes(state, values, fill=state)) + ggplot2::geom_boxplot()+
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free')

