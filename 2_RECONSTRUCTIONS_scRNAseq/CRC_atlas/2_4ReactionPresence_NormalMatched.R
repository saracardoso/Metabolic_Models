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

load('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/general/reports/NormalMatched_env.RData')






# --------------------------------------------------------
# --- Percentage of active reactions in the subsystems ---
# --------------------------------------------------------

# 0. Load annotations (only considering pathways with more than 5 reactions):
HumanGEM_subsystems_rxns = jsonlite::read_json('./GENERAL/utility_data/subsystems_reactions_mapping.json',
                                               simplifyVector = TRUE)
for(path in names(HumanGEM_subsystems_rxns)){
  rxns = HumanGEM_subsystems_rxns[[path]][HumanGEM_subsystems_rxns[[path]]%in%rownames(rxn_presence_gprs)]
  if(length(rxns)>5) HumanGEM_subsystems_rxns[[path]] = rxns
  else HumanGEM_subsystems_rxns[[path]] = NULL
}

# 1. Calculate percentage of active reactions per pathway per model:
n_models = dim(rxn_presence)[2]
n_pathways = length(HumanGEM_subsystems_rxns)
# 1.1. Create matrix where results are stored:
perc_presence = matrix(rep(0, n_models*n_pathways), nrow=n_pathways)
dimnames(perc_presence) = list(names(HumanGEM_subsystems_rxns), colnames(rxn_presence))
# 1.2. Populate matrix:
for(pathway in rownames(perc_presence)){
  message(pathway)
  for(model in colnames(perc_presence)){
    presence_pathway_model = rxn_presence[HumanGEM_subsystems_rxns[[pathway]],model]
    perc = sum(presence_pathway_model) / length(presence_pathway_model) * 100
    perc_presence[pathway, model] = perc
  }
}
# 1.3. Heatmap:
ht = ComplexHeatmap::Heatmap(perc_presence, col = circlize::colorRamp2(c(0, 50, 100), c('blue', 'white', 'red')),
                        name='% Presence',
                        column_split = metadata[colnames(perc_presence),]$cell_type, column_title = NULL,
                        show_column_names = FALSE, row_names_gp = grid::gpar(fontsize=6),
                        top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata[colnames(perc_presence),]$cell_type,
                                                                           state=metadata[colnames(perc_presence),]$state,
                                                                           col=list(cell_type=ct_colors, state=state_colors),
                                                                           annotation_legend_param = list(cell_type=list(title='Cell Type'),
                                                                                                          state=list(title='State'))))
ComplexHeatmap::draw(ht, merge_legend = TRUE)
# 1.4. Differential pathways:
metadata = metadata[colnames(perc_presence),]
diff_pathways = c()
for(ct in unique(metadata$cell_type)){
  message(ct)
  ct_paths_means = rowMeans(perc_presence[,metadata$cell_type==ct])
  other_paths_means = rowMeans(perc_presence[,metadata$cell_type!=ct])
  fc = ct_paths_means / other_paths_means
  fc[fc<1] = (-1/fc[fc<1])
  names(fc) = row.names(perc_presence)
  paths_test = names(fc)[fc>1 | fc<(-1)]
  test_meta1 = metadata$cell_type
  test_meta1 = data.frame(cell_type=test_meta1, row.names=colnames(perc_presence))
  test_meta1$cell_type[test_meta1$cell_type!=ct] = 'other'
  p_val = sapply(paths_test,
                 FUN = function(x){
                   return(wilcox.test(as.numeric(perc_presence[x,]) ~ as.vector(test_meta1[,'cell_type']))$p.value)
                 })
  pval_adjust = p.adjust(p_val, method='fdr')
  res_ct = data.frame(pathway=paths_test, pval_adjust, pval=p_val, FC=fc[paths_test],cell_type=ct,
                      row.names=paste(paths_test, ct, sep='_'))
  diff_pathways = rbind(diff_pathways, res_ct)
}
write.csv(diff_pathways, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/diff_analysis_percPresence.csv')
# 1.5. Heatmap with only differential pathways:
diff_paths_heatmap = unique(diff_pathways$pathway[diff_pathways$pval_adjust<0.01])
ht2 = ComplexHeatmap::Heatmap(perc_presence[diff_paths_heatmap,], col = circlize::colorRamp2(c(0, 50, 100), c('blue', 'white', 'red')),
                             name='% Presence',
                             column_split = metadata[colnames(perc_presence),]$cell_type, column_title = NULL,
                             show_column_names = FALSE, row_names_gp = grid::gpar(fontsize=6),
                             top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata[colnames(perc_presence),]$cell_type,
                                                                                state=metadata[colnames(perc_presence),]$state,
                                                                                col=list(cell_type=ct_colors, state=state_colors),
                                                                                annotation_legend_param = list(cell_type=list(title='Cell Type'),
                                                                                                               state=list(title='State'))))
ComplexHeatmap::draw(ht2, merge_legend = TRUE)
