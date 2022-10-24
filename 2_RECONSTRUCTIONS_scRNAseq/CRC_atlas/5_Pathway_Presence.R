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

cms_colors = c('#d4cfcf', '#cc79a7', '#b3c9e6', '#c3d7a4', '#52854c', '#fdb981')
names(cms_colors) = c('CMS3', 'Mixed', 'CMS1', 'CMS2', 'Normal Matched', 'CMS4')

# Load data:
rxn_presence = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/Genes_to_Fluxes/reaction_presence.csv',
                             row.names=1, check.names=FALSE)
metadata = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/metadata.csv',
                    row.names=1, check.names=FALSE)






# --------------------------------------------------------
# --- Percentage of active reactions in the subsystems ---
# --------------------------------------------------------


# 0. Load annotations (only considering pathways with more than 5 reactions):
HumanGEM_subsystems_rxns = jsonlite::read_json('./GENERAL/utility_data/subsystems_reactions_mapping.json',
                                               simplifyVector = TRUE)
for(path in names(HumanGEM_subsystems_rxns)){
  rxns = HumanGEM_subsystems_rxns[[path]][HumanGEM_subsystems_rxns[[path]]%in%rownames(rxn_presence)]
  if(length(rxns)>5) HumanGEM_subsystems_rxns[[path]] = rxns
  else HumanGEM_subsystems_rxns[[path]] = NULL
}


# 1. Calculate percentage of active reactions per pathway per model:
n_models = dim(metadata)[1]
n_pathways = length(HumanGEM_subsystems_rxns)

# 1.1. Create matrix where results are stored:
perc_presence = matrix(rep(0, n_models*n_pathways), nrow=n_pathways)
dimnames(perc_presence) = list(names(HumanGEM_subsystems_rxns), rownames(metadata))

# 1.2. Populate matrix:
for(pathway in rownames(perc_presence)){
  message(pathway)
  for(model in colnames(perc_presence)){
    presence_pathway_model = rxn_presence[HumanGEM_subsystems_rxns[[pathway]],model]
    perc = sum(presence_pathway_model) / length(presence_pathway_model) * 100
    perc_presence[pathway, model] = perc
  }
}

# 1.3. Heatmap with all pathways:
ht = ComplexHeatmap::Heatmap(perc_presence, col = circlize::colorRamp2(c(0, 50, 100), c('blue', 'white', 'red')),
                        name='% Presence',
                        column_split = metadata[colnames(perc_presence),]$cell_type, column_title = NULL,
                        show_column_names = FALSE, row_names_gp = grid::gpar(fontsize=6),
                        top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata[colnames(perc_presence),]$cell_type,
                                                                           state=metadata[colnames(perc_presence),]$state,
                                                                           cms=metadata[colnames(perc_presence),]$CMS,
                                                                           col=list(cell_type=ct_colors, state=state_colors,
                                                                                    cms=cms_colors),
                                                                           annotation_legend_param = list(cell_type=list(title='Cell Type'),
                                                                                                          state=list(title='State'),
                                                                                                          cms=list(title='CMS'))))
ComplexHeatmap::draw(ht, merge_legend = TRUE)

# 1.4. Pathways present/not present across most cell-types:
paths_present = c()
paths_not_present = c()
for(path in rownames(perc_presence)){
  if(summary(perc_presence[path,])[3] > 75 & (sum(perc_presence[path,] > 75)/dim(perc_presence)[2] > .80)) paths_present = c(paths_present, path)
  if(summary(perc_presence[path,])[3] < 25 & (sum(perc_presence[path,] < 25)/dim(perc_presence)[2] > .80)) paths_not_present = c(paths_not_present, path)
}
ht = ComplexHeatmap::Heatmap(perc_presence[c(paths_present, paths_not_present), ],
                             col = circlize::colorRamp2(c(0, 50, 100), c('blue', 'white', 'red')),
                             name='% Presence',
                             #row_names_rot = 325,
                             row_names_max_width = ComplexHeatmap::max_text_width(rownames(perc_presence[c(paths_present, paths_not_present), ]), 
                                                                  gp = grid::gpar(fontsize = 12)),
                             row_split = c(rep('present', length(paths_present)), rep('not_present', length(paths_not_present))), row_title = NULL,
                             column_split = metadata[colnames(perc_presence),]$cell_type, column_title = NULL,
                             show_column_names = FALSE, row_names_gp = grid::gpar(fontsize=8),
                             top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata[colnames(perc_presence),]$cell_type,
                                                                                state=metadata[colnames(perc_presence),]$state,
                                                                                cms=metadata[colnames(perc_presence),]$CMS,
                                                                                col=list(cell_type=ct_colors, state=state_colors,
                                                                                         cms=cms_colors),
                                                                                annotation_legend_param = list(cell_type=list(title='Cell Type'),
                                                                                                               state=list(title='State'),
                                                                                                               cms=list(title='CMS')),
                                                                                annotation_label = c('Cell Type', 'State', 'CMS')))
ComplexHeatmap::draw(ht, merge_legend = TRUE, heatmap_legend_side = 'bottom', show_heatmap_legend=FALSE)

# 1.5. Pathways differentially present between normal and tumour for each cell-type:
diff_pathways = c()
for(ct in c("Cytotoxic CD8 Tcells", "IL17+ CD4 Tcells", "Memory CD4 Tcells", "Memory CD8 Tcells", "Naive CD4 Tcells", "Naive CD8 Tcells",
            "Regulatory CD4 Tcells")){
  message(ct)
  groups_vector = metadata$state[metadata$cell_type==ct]
  cts_vector = rownames(metadata)[metadata$cell_type==ct]
  
  norm_cts = cts_vector[groups_vector=='Normal Matched']
  tum_cts = cts_vector[groups_vector!='Normal Matched']
  paths_to_test = c()
  fcs_paths = c()
  for(path in rownames(perc_presence)){
    fc = gtools::foldchange(mean(perc_presence[path, tum_cts]), mean(perc_presence[path, norm_cts]))
    if(!is.na(fc) & abs(fc) > 1.5){
      paths_to_test = c(paths_to_test, path)
      fcs_paths = c(fcs_paths, fc)
    }
  }
  if(length(paths_to_test) != 0){
    p_val = sapply(paths_to_test,
                   FUN = function(x){
                     return(wilcox.test(as.numeric(perc_presence[x, cts_vector]) ~ groups_vector)$p.value)
                   })
    pval_adjust = p.adjust(p_val, method='fdr')
    res_ct = data.frame(pathway=paths_to_test, pval_adjust, pval=p_val, fold_change=fcs_paths, cell_type=ct)
    diff_pathways = rbind(diff_pathways, res_ct)
  }
}
write.csv(diff_pathways, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/pathways_presence/normalVStumour.csv')

# 1.6. Heatmap for each cell-type of the above results:
hetamaps_normVStum = list()
for(ct in unique(diff_pathways$cell_type)){
  ct_diffs = na.omit(diff_pathways$pathway[diff_pathways$pval_adjust < 0.05 & diff_pathways$cell_type==ct])
  ct_models = rownames(metadata)[metadata$cell_type==ct]
  if(length(ct_diffs) == 1){
    mat = t(as.data.frame(perc_presence[ct_diffs, ct_models]))
    rownames(mat) = ct_diffs
  } 
  else mat = perc_presence[ct_diffs, ct_models]
  hetamaps_normVStum[[ct]] = ComplexHeatmap::Heatmap(mat,
                               col = circlize::colorRamp2(c(0, 50, 100), c('blue', 'white', 'red')),
                               name='% Presence',
                               #row_names_rot = 325,
                               column_km = 3,
                               row_names_max_width = ComplexHeatmap::max_text_width(rownames(perc_presence[ct_diffs, ]), 
                                                                                    gp = grid::gpar(fontsize = 12)),
                               show_column_names = FALSE, row_names_gp = grid::gpar(fontsize=8),
                               top_annotation = ComplexHeatmap::HeatmapAnnotation(state=metadata[ct_models,]$state,
                                                                                  cms=metadata[ct_models,]$CMS,
                                                                                  col=list(state=state_colors, cms=cms_colors),
                                                                                  annotation_legend_param = list(state=list(title='State'),
                                                                                                                 cms=list(title='CMS')),
                                                                                  annotation_label = c('State', 'CMS')),
                               heatmap_legend_param = list(direction = "horizontal"))
}
ComplexHeatmap::draw(hetamaps_normVStum$`Cytotoxic CD8 Tcells`, merge_legend = TRUE, heatmap_legend_side = 'bottom', show_heatmap_legend=TRUE)
ComplexHeatmap::draw(hetamaps_normVStum$`Regulatory CD4 Tcells`, merge_legend = TRUE, heatmap_legend_side = 'bottom', show_heatmap_legend=FALSE)
ComplexHeatmap::draw(hetamaps_normVStum$`IL17+ CD4 Tcells`, merge_legend = TRUE, heatmap_legend_side = 'bottom', show_heatmap_legend=TRUE)

# 1.7. Pathways differentially present between pairs of cell-types:
diff_pathways = c()
ct_pairs = t(combn(names(ct_colors), 2))
for(idx in 1:dim(ct_pairs)[1]){
  ct_pair = ct_pairs[idx,]
  message(paste(idx, ct_pair, collapse='_'))
  models = rownames(metadata)[metadata$cell_type%in%ct_pair]
  groups_vector = metadata[models, 'cell_type']
  
  ct1 = rownames(metadata)[metadata$cell_type==ct_pair[1]]
  ct2 = rownames(metadata)[metadata$cell_type==ct_pair[2]]
  paths_to_test = c()
  fcs_paths = c()
  for(path in rownames(perc_presence)){
    fc = gtools::foldchange(mean(perc_presence[path, ct1]), mean(perc_presence[path, ct2]))
    if(!is.na(fc) & abs(fc) > 2){
      paths_to_test = c(paths_to_test, path)
      fcs_paths = c(fcs_paths, fc)
    }
  }
  if(length(paths_to_test) != 0){
    p_val = sapply(paths_to_test,
                   FUN = function(x){
                     return(wilcox.test(as.numeric(perc_presence[x, models]) ~ groups_vector)$p.value)
                   })
    pval_adjust = p.adjust(p_val, method='fdr')
    res_ct = data.frame(pathway=paths_to_test, pval_adjust, pval=p_val, fold_change=fcs_paths, cell_type1=ct_pair[1],
                        cell_type2=ct_pair[2])
    diff_pathways = rbind(diff_pathways, res_ct)
  }
}
write.csv(diff_pathways, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/pathways_presence/ct_diffPathways.csv')

# 1.6. Heatmap of the above results:
heatmaps_ctpairs = list()
for(idx in 1:dim(ct_pairs)[1]){
  ct_pairs[idx,]
  ct_diffs_paths = diff_pathways[diff_pathways$cell_type1 == ct_pairs[idx, 1] & diff_pathways$cell_type2 == ct_pairs[idx, 2] &
                                   diff_pathways$pval_adjust < 0.05,]
  ct_diffs_models = rownames(metadata)[metadata$cell_type%in%ct_pairs[idx,]]
  if(length(ct_diffs_paths$pathway) == 1){
    mat = t(as.data.frame(perc_presence[ct_diffs_paths$pathway, ct_diffs_models]))
    rownames(mat) = ct_diffs_paths$pathway
  } 
  else mat = perc_presence[ct_diffs_paths$pathway, ct_diffs_models]
  heatmaps_ctpairs[[idx]] = ComplexHeatmap::Heatmap(mat,
                               col = circlize::colorRamp2(c(0, 50, 100), c('blue', 'white', 'red')),
                               name='% Presence',
                               #row_names_rot = 325,
                               row_names_max_width = ComplexHeatmap::max_text_width(rownames(mat), 
                                                                                    gp = grid::gpar(fontsize = 12)),
                               #column_split = metadata[ct_diffs_models,]$cell_type, column_title = NULL,
                               show_column_names = FALSE, row_names_gp = grid::gpar(fontsize=8),
                               top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata[ct_diffs_models,]$cell_type,
                                                                                  state=metadata[ct_diffs_models,]$state,
                                                                                  cms=metadata[ct_diffs_models,]$CMS,
                                                                                  col=list(cell_type=ct_colors, state=state_colors, cms=cms_colors),
                                                                                  annotation_legend_param = list(cell_type=list(title='Cell Type'),
                                                                                                                 state=list(title='State'),
                                                                                                                 cms=list(title='CMS')),
                                                                                  annotation_label = c('Cell Type', 'State', 'CMS')),
                               heatmap_legend_param = list(direction = "horizontal"))
}
# 1.6.1. Focus on the pairs:
ComplexHeatmap::draw(heatmaps_ctpairs[[24]], merge_legend = TRUE, heatmap_legend_side = 'bottom')#, show_heatmap_legend=FALSE)
ComplexHeatmap::draw(heatmaps_ctpairs[[36]], merge_legend = TRUE, heatmap_legend_side = 'bottom')#, show_heatmap_legend=FALSE)
ComplexHeatmap::draw(heatmaps_ctpairs[[41]], merge_legend = TRUE, heatmap_legend_side = 'bottom')#, show_heatmap_legend=FALSE)
# 1.6.2. Putting all significant pathways together for all cell-types:
paths_plot = unique(diff_pathways[diff_pathways$pval_adjust < 0.05, 'pathway'])
mat = perc_presence[paths_plot, rownames(metadata)]
ht = ComplexHeatmap::Heatmap(mat,
                             col = circlize::colorRamp2(c(0, 50, 100), c('blue', 'white', 'red')),
                             name='% Presence',
                             #row_names_rot = 325,
                             row_names_max_width = ComplexHeatmap::max_text_width(rownames(mat), 
                                                                                  gp = grid::gpar(fontsize = 12)),
                             column_split = metadata$cell_type, column_title = NULL,
                             show_column_names = FALSE, row_names_gp = grid::gpar(fontsize=8),
                             top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata$cell_type,
                                                                                state=metadata$state,
                                                                                cms=metadata$CMS,
                                                                                col=list(cell_type=ct_colors, state=state_colors, cms=cms_colors),
                                                                                annotation_legend_param = list(cell_type=list(title='Cell Type'),
                                                                                                               state=list(title='State'),
                                                                                                               cms=list(title='CMS')),
                                                                                annotation_label = c('Cell Type', 'State', 'CMS')),
                             heatmap_legend_param = list(direction = "horizontal"))
ComplexHeatmap::draw(ht, merge_legend = TRUE, heatmap_legend_side = 'bottom')#, show_heatmap_legend=FALSE)






# ----------------------------
# --- Structure similarity ---
# ----------------------------

distance_matrix = as.matrix(dist(t(rxn_presence)))

ht = ComplexHeatmap::Heatmap(distance_matrix,
                             col = circlize::colorRamp2(c(0, 35, 70), c('blue', 'white', 'red')),
                             name='Euclidean Distance',
                             column_split = metadata$cell_type, column_title = NULL,
                             row_split = metadata$cell_type, row_title = NULL,
                             show_column_names = FALSE, show_row_names = FALSE,
                             top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata$cell_type,
                                                                                state=metadata$state,
                                                                                cms=metadata$CMS,
                                                                                col=list(cell_type=ct_colors, state=state_colors, cms=cms_colors),
                                                                                annotation_legend_param = list(cell_type=list(title='Cell Type'),
                                                                                                               state=list(title='State'),
                                                                                                               cms=list(title='CMS')),
                                                                                annotation_label = c('Cell Type', 'State', 'CMS')),
                             left_annotation = ComplexHeatmap::rowAnnotation(cell_type=metadata$cell_type,
                                                                             state=metadata$state,
                                                                             cms=metadata$CMS,
                                                                             col=list(cell_type=ct_colors, state=state_colors, cms=cms_colors),
                                                                             annotation_legend_param = list(cell_type=list(title='Cell Type'),
                                                                                                            state=list(title='State'),
                                                                                                            cms=list(title='CMS')),
                                                                             annotation_label = c('Cell Type', 'State', 'CMS')),
                             heatmap_legend_param = list(direction = "horizontal"))
ComplexHeatmap::draw(ht, merge_legend = TRUE, heatmap_legend_side = 'right')#, show_heatmap_legend=FALSE)
