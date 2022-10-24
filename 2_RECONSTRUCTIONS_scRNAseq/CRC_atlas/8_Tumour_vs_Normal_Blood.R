ct_colors = c('#b30000', '#999999', '#cc79a7', '#b3c9e6', '#4e84c4', '#c3d7a4', '#52854c', '#fdb981',
              '#d16103', '#293352')
names(ct_colors) = c('Cytotoxic CD8 Tcells', 'Follicular CD4 Tcells', 'IL17+ CD4 Tcells',
                     'Memory CD4 Tcells', 'Memory CD8 Tcells', 'Naive CD4 Tcells', 'Naive CD8 Tcells',
                     'Proliferative CD4 Tcells', 'Proliferative CD8 Tcells',  'Regulatory CD4 Tcells')
state_colors = c('#52854c', '#b30000', '#ff6666', '#cc79a7')
names(state_colors) = c('Normal Matched', 'Tumour', 'Tumour Core', 'Tumour Border')



# -----------------
# --- Read Data ---
# -----------------

normal_medium = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/Normal_Blood_pFBA.csv',
                         row.names=1, header=TRUE, check.names=FALSE)
normal_medium = 2 * ((1 / (1 + exp(-normal_medium))) - 0.5)
tumour_medium = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/Tumour_Blood_pFBA.csv',
                         row.names=1, header=TRUE, check.names=FALSE)
tumour_medium = 2 * ((1 / (1 + exp(-tumour_medium))) - 0.5)

# Metadata:
metadata = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/metadata.csv',
                    row.names=1, check.names=FALSE)





# --------------------------------------------------------------------------
# --- Do models' fluxes change between using a normal and tumour medium? ---
# --------------------------------------------------------------------------

models_analysed = intersect(colnames(normal_medium), colnames(tumour_medium))

# 1. Perform the statistical tests:
test_results_nvst = data.frame(p_value=rep(NA, length(models_analysed)),
                               p_adjust=rep(NA, length(models_analysed)))
row.names(test_results_nvst) = models_analysed#colnames(tumour_medium)
pb <- txtProgressBar(min=0, max=length(models_analysed), style=3)
i=0
for(model in models_analysed){#colnames(tumour_medium)){
  i = i+1
  to_test = data.frame(medium=c(rep('Normal', dim(tumour_medium)[1]),
                                rep('Tumour', dim(tumour_medium)[1])),
                       flux=c(as.numeric(normal_medium[, model]),
                              as.numeric(tumour_medium[, model])))
  wilcox_test_result = wilcox.test(to_test$flux ~ to_test$medium, paired=TRUE)
  test_results_nvst[model,'p_value'] = wilcox_test_result$p.value
  setTxtProgressBar(pb, i)
}
close(pb)
test_results_nvst$p_adjust = p.adjust(test_results_nvst$p_value, method='fdr')
test_results_nvst = na.omit(test_results_nvst)
write.csv(test_results_nvst, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/TumourVSNormal/normal_vs_tumour_pFBA.csv')

# 2. Percentage of models that differ
sum(test_results_nvst$p_adjust < 0.05) / dim(test_results_nvst)[1] * 100 # 53 models differ significantly, i.e., ~46%

# 3. Percentage of models that differ, per cell-type
differ_per_ct = table(metadata[rownames(test_results_nvst), 'cell_type'][test_results_nvst$p_adjust < 0.05])[unique(metadata$cell_type)]
differ_per_ct[is.na(differ_per_ct)] = 0
names(differ_per_ct) = unique(metadata$cell_type)
x =  differ_per_ct / table(metadata[rownames(test_results_nvst), 'cell_type'])[unique(metadata$cell_type)] * 100

# 4. Plot this information
perc_df = data.frame(fdr=test_results_nvst$p_adjust,
                     difference=test_results_nvst$p_adjust < 0.05,
                     cell_type = metadata[rownames(test_results_nvst), 'cell_type'],
                     state = metadata[rownames(test_results_nvst), 'state'])
perc_df$state[grep('Tumour', perc_df$state)] = 'Tumour'

perc_df$cell_type = factor(perc_df$cell_type, levels=names(ct_colors))

perc_df2 = data.frame(difference=as.numeric(x), cell_type = names(x))

ggplot2::ggplot(perc_df[perc_df$difference,], ggplot2::aes(cell_type, fill=cell_type)) +
  ggplot2::geom_bar() + ggplot2::scale_fill_manual(values=ct_colors) +
  ggplot2::theme(axis.text.x=ggplot2::element_blank(), legend.position='none') +
  #ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=7)) +
  ggplot2::ylab('Number of Models') + ggplot2::xlab('Cell Types')

ggplot2::ggplot(perc_df2, ggplot2::aes(x=cell_type, y=difference, fill=cell_type)) +
  ggplot2::geom_bar(stat='identity') + ggplot2::scale_fill_manual(values=ct_colors) +
  ggplot2::theme(legend.position='none', axis.text.x=ggplot2::element_blank()) +
  ggplot2::ylab('Percentage') + ggplot2::xlab('Cell Types') + ggplot2::ylim(0, 80)

ggplot2::ggplot(perc_df, ggplot2::aes(y=fdr, x=cell_type)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::scale_colour_manual(values=state_colors[c('Normal Matched', 'Tumour')]) +
  ggplot2::geom_jitter(ggplot2::aes(colour=state)) + ggplot2::geom_hline(yintercept=0.05, colour='grey') +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=7)) +
  ggplot2::ylab('FDR') + ggplot2::xlab('Cell Types')





# -------------------------------------------------------------------
# --- How much subsystems differ between normal and tumour medium ---
# -------------------------------------------------------------------

# 0. Models that significantly differ:
diff_models = rownames(test_results_nvst)[test_results_nvst$p_adjust < 0.05]

# 1. For each model that differs, evaluate the differences:
reaction_differences = tumour_medium[, diff_models] - normal_medium[, diff_models]
write.csv(reaction_differences,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/TumourVSNormal/rxn_differences_pFBA.csv')
reaction_FC = tumour_medium[, diff_models] / normal_medium[, diff_models]
write.csv(reaction_FC,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/TumourVSNormal/rxn_FC_pFBA.csv')

# 2. Metadata of only the differing models:
metadata_diff = metadata[diff_models, ]
metadata_diff$if_diff = rep('No', length(diff_models))
metadata_diff[diff_models, 'if_diff'] = 'Yes'

# 3. For each model, check what reactions differ:
# 3.1. Get mapping between reactions and pathways:
rxns_subsystems = jsonlite::read_json('./GENERAL/utility_data/reactions_subsystems_mapping.json',
                                      simplifyVector=TRUE)
rxns = c()
paths = c()
for(rxn in rownames(reaction_differences)){
  rxns = c(rxns, rep(rxn, length(rxns_subsystems[[rxn]])))
  paths = c(paths, rxns_subsystems[[rxn]])
}
# 3.2. Get how reactions change from normal to tumour:
rxn_changes = list()
i=0
for(model in rownames(metadata_diff)){
  i=i+1
  message(model, ' | ', i)
  rxn_changes[[model]] = data.frame(reaction = rxns, state=rep('Unchanged', length(rxns)), pathway=paths)
  for(idx in 1:length(rxns)){
    rxn = rxns[idx]
    if((reaction_differences[rxn, model] != 0)){ #& 
        #(abs(reaction_FC[rxn, model]) > 1.5 | abs(1/reaction_FC[rxn, model]) > 1.5)){
      if(reaction_FC[rxn, model] < 0) rxn_changes[[model]][idx, 'state'] = 'Reversed'
      else{
        if(abs(tumour_medium[rxn, model]) > abs(normal_medium[rxn, model]))
          rxn_changes[[model]][idx, 'state'] = 'Increased'
        else rxn_changes[[model]][idx, 'state'] = 'Decreased'
      }
    }
  }
}
jsonlite::write_json(rxn_changes,
                     './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/TumourVSNormal/rxn_changes_pFBA.json')
# 3.3. Get ratio of reactions in subsystems that changed:
# 3.3.1. Get the percentages:
perc_changed_paths = matrix(rep(0, length(unique(paths))*length(rxn_changes)),
                            ncol=length(rxn_changes))
rownames(perc_changed_paths) = unique(paths)
colnames(perc_changed_paths) = names(rxn_changes)
for(model in names(rxn_changes)){
  for(path in unique(paths)){
    rxns_path = rxn_changes[[model]]$pathway == path
    n_rxns_changed = sum(rxn_changes[[model]]$state[rxns_path] != 'Unchanged')
    perc_changed_paths[path, model] = n_rxns_changed / sum(rxns_path)
  }
}
# 3.3.2. Median percentage per pathway across all models:
median_perc_models = c()
for(path in rownames(perc_changed_paths)){
  median_perc_models = c(median_perc_models, median(perc_changed_paths[path, ]))
}
names(median_perc_models) = rownames(perc_changed_paths)
top30_most = sort(median_perc_models, decreasing=T)[1:30] # Top 30 most affected pathways
# 3.4. Boxplot of top30 most affected pathways:
percs = c()
paths = c()
mods = c()
cts = c()
for(model in colnames(perc_changed_paths)){
  percs = c(percs, perc_changed_paths[names(top30_most), model])
  paths = c(paths, names(top30_most))
  mods = c(mods, rep(model, 30))
  cts = c(cts, rep(metadata$cell_type[rownames(metadata)==model], 30))
}
df_paths = data.frame(percs, paths, mods, cts)
df_paths$paths = factor(df_paths$paths, levels=names(top30_most))

ggplot2::ggplot(df_paths, ggplot2::aes(x=percs, y=paths)) +
  ggplot2::geom_boxplot() +
  ggplot2::xlab('Ratio Affected Reactions') + ggplot2::ylab('Pathways')

# 3.5. Visualize the percentages of these pathways:
mat = perc_changed_paths[c(names(top30_most)), rownames(metadata_diff)]
heat_perc = ComplexHeatmap::Heatmap(mat,
                                    row_names_max_width = ComplexHeatmap::max_text_width(rownames(mat), 
                                                                                         gp = grid::gpar(fontsize = 12)),
                                  column_split = metadata_diff$cell_type,
                                  #row_split = c(rep('Least', 15), rep('Most', 15)),
                                  column_title = NULL, row_title = NULL,
                                  show_column_names = FALSE,
                                  row_names_gp = grid::gpar(fontsize=7),
                                  heatmap_legend_param = list(title='Ratio Affected reactions',
                                                              title_position = "lefttop-rot"),
                                  cluster_rows = FALSE,
                                  top_annotation = ComplexHeatmap::HeatmapAnnotation(
                                    cell_type=metadata_diff$cell_type,
                                    state=metadata_diff$state,
                                    if_diff=metadata_diff$if_diff,
                                    col=list(cell_type=ct_colors, state=state_colors,
                                             if_diff=c('Yes'='lightgreen', 'No'='grey')),
                                    annotation_legend_param = list(cell_type=list(title='Cell Type',
                                                                                  title_position = "lefttop-rot"),
                                                                   state=list(title='State',
                                                                              title_position = "lefttop-rot"),
                                                                   if_diff=list(title='N-T Difference',
                                                                                title_position = "lefttop-rot")),
                                    show_annotation_name=FALSE, simple_anno_size = grid::unit(2,'mm')))
ComplexHeatmap::draw(heat_perc, merge_legend = TRUE, heatmap_legend_side='bottom',
                     annotation_legend_side='bottom')

# 3.6. How biomass changes between normal and tumour medium:
#plot_df = data.frame(flux=c(as.numeric(normal_medium['MAR13082', rownames(metadata)]),
#                            as.numeric(tumour_medium['MAR13082', rownames(metadata)])),
#                     medium=c(rep('Normal', dim(metadata)[1]), rep('Tumour', dim(metadata)[1])),
#                     ct=metadata$cell_type, state=metadata$state, CMS=metadata$CMS,
#                     paired = c(1:dim(metadata)[1], 1:dim(metadata)[1]))
#ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=flux, colour=ct)) +
#  ggplot2::geom_boxplot(outlier.shape=NA) +
#  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey', position = ggplot2::position_dodge(0.15)) +
#  ggplot2::facet_wrap(ggplot2::vars(ct)) + 
#  #ggplot2::geom_jitter(width = 0.15) +
#  ggplot2::geom_point(ggplot2::aes(colour=ct, group=paired),
#                      position = ggplot2::position_dodge(0.2)) +
#  ggplot2::theme_light() + 
#  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8), legend.position='none') +
#  ggplot2::scale_colour_manual(values=ct_colors) +
#  ggplot2::xlab('') + ggplot2::ylab('Flux (mmol/gDW/h)')



# 3.6. How biomass changes between normal and tumour medium, only for differing models:
models_diff = rownames(test_results_nvst)[test_results_nvst$p_adjust < 0.05]
plot_df = data.frame(flux=c(as.numeric(normal_medium['MAR13082', models_diff]),
                            as.numeric(tumour_medium['MAR13082', models_diff])),
                     medium=c(rep('Normal', length(models_diff)), rep('Tumour', length(models_diff))),
                     ct=metadata[models_diff, 'cell_type'], state=metadata[models_diff, 'state'], CMS=metadata[models_diff, 'CMS'],
                     paired = c(1:length(models_diff), 1:length(models_diff)))
ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=flux, colour=ct)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey', position = ggplot2::position_dodge(0.15)) +
  ggplot2::facet_wrap(ggplot2::vars(ct)) + 
  #ggplot2::geom_jitter(width = 0.15) +
  ggplot2::geom_point(ggplot2::aes(colour=ct, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') + ggplot2::ylab('Flux (mmol/gDW/h)')


# 3.7. How biomass changed
biomass_change = rep('No Change', length(models_analysed))
names(biomass_change) = models_analysed

diff_all = normal_medium['MAR13082', models_analysed] - tumour_medium['MAR13082', models_analysed]
sum(diff_all != 0)

models_diff_biomass = names(diff_all)[diff_all > 0]
#sort(table(metadata[models_diff_biomass, 'cell_type']), decreasing = T)
biomass_change[models_diff_biomass] = 'Decrease'

models_diff_biomass = names(diff_all)[diff_all < 0]
#sort(table(metadata[models_diff_biomass, 'cell_type']), decreasing = T)
biomass_change[models_diff_biomass] = 'Increase'

colors_changes = c('#b3c9e6', '#b30000', '#52854c')
names(colors_changes) = c('No Change', 'Increase', 'Decrease')
plot_df = data.frame(biomass_change, cell_type=metadata[models_analysed, 'cell_type'],
                     state=metadata[models_analysed, 'state'],
                     cms=metadata[models_analysed, 'CMS'])
plot_df$biomass_change = factor(plot_df$biomass_change, levels=c('No Change', 'Increase', 'Decrease'))
ggplot2::ggplot(plot_df, ggplot2::aes(biomass_change, fill=biomass_change)) +
  ggplot2::geom_bar(position='stack') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8), legend.position='none')+
  ggplot2::scale_fill_manual(values=colors_changes) +
  ggplot2::xlab('') + ggplot2::ylab('')
ggplot2::ggplot(plot_df, ggplot2::aes(cell_type, fill=biomass_change)) +
  ggplot2::geom_bar(position='stack') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8), legend.position='none')+
  ggplot2::scale_fill_manual(values=colors_changes) +
  ggplot2::xlab('') + ggplot2::ylab('')
ggplot2::ggplot(plot_df, ggplot2::aes(state, fill=biomass_change)) +
  ggplot2::geom_bar(position='stack') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8), legend.position='none')+
  ggplot2::scale_fill_manual(values=colors_changes) +
  ggplot2::xlab('') + ggplot2::ylab('')
ggplot2::ggplot(plot_df, ggplot2::aes(cms, fill=biomass_change)) +
  ggplot2::geom_bar(position='stack') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_fill_manual(values=colors_changes) +
  ggplot2::xlab('') + ggplot2::ylab('')



# 3.4. Distribution of the types of change in the top 15 subsystems:
df_plot = c()
for(model in names(rxn_changes)){
  model_df = rxn_changes[[model]]
  temp_df = model_df[model_df$state!='Unchanged' & model_df$pathway%in%names(top30_most),
                     c('reaction', 'state', 'pathway')]
  temp_df$cell_type = metadata_diff[model, 'cell_type']
  temp_df$model = model
  if(!is.null(df_plot)) df_plot = rbind(df_plot, temp_df)
  else df_plot = temp_df
}
ggplot2::ggplot(df_plot[df_plot$cell_type=='IL17+ CD4 Tcells',], ggplot2::aes(model, fill=state)) + 
  ggplot2::geom_bar(position="fill") +
  ggplot2::scale_fill_manual(values=c('Increased'='darkgreen', 'Decreased'='darkred', 'Reversed'='grey')) +
  ggplot2::theme(axis.text.x=ggplot2::element_blank()) +
  ggplot2::facet_wrap(ggplot2::vars(pathway), ncol=3)

ggplot2::ggplot(df_plot[df_plot$cell_type=='Proliferative CD8 Tcells',], ggplot2::aes(model, fill=state)) + 
  ggplot2::geom_bar(position="fill") +
  ggplot2::scale_fill_manual(values=c('Increased'='darkgreen', 'Decreased'='darkred', 'Reversed'='grey')) +
  ggplot2::theme(axis.text.x=ggplot2::element_blank()) +
  ggplot2::facet_wrap(ggplot2::vars(pathway), ncol=3)


ggplot2::ggplot(df_plot[df_plot$pathway=='Glycolysis / Gluconeogenesis',], ggplot2::aes(model, fill=state)) + 
  ggplot2::geom_bar(position="fill") +
  ggplot2::scale_fill_manual(values=c('Increased'='darkgreen', 'Decreased'='darkred', 'Reversed'='grey')) +
  ggplot2::theme(axis.text.x=ggplot2::element_blank()) +
  ggplot2::facet_wrap(ggplot2::vars(cell_type), scales='free')

ggplot2::ggplot(df_plot[df_plot$pathway=='Oxidative phosphorylation',], ggplot2::aes(model, fill=state)) + 
  ggplot2::geom_bar(position="fill") +
  ggplot2::scale_fill_manual(values=c('Increased'='darkgreen', 'Decreased'='darkred', 'Reversed'='grey')) +
  ggplot2::theme(axis.text.x=ggplot2::element_blank()) +
  ggplot2::facet_wrap(ggplot2::vars(cell_type), scales='free')

ggplot2::ggplot(df_plot[df_plot$pathway=='Beta oxidation of unsaturated fatty acids (n-7) (mitochondrial)',], ggplot2::aes(model, fill=state)) + 
  ggplot2::geom_bar(position="fill") +
  ggplot2::scale_fill_manual(values=c('Increased'='darkgreen', 'Decreased'='darkred', 'Reversed'='grey')) +
  ggplot2::theme(axis.text.x=ggplot2::element_blank()) +
  ggplot2::facet_wrap(ggplot2::vars(cell_type), scales='free')



