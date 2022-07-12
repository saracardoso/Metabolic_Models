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

normal_medium = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/normal_FBA.csv',
                         row.names=1, header=TRUE, check.names=FALSE)
normal_medium = 2 * ((1 / (1 + exp(-normal_medium))) - 0.5)
tumour_medium = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/4_tumourMedium/normal_FBA.csv',
                         row.names=1, header=TRUE, check.names=FALSE)
tumour_medium = 2 * ((1 / (1 + exp(-tumour_medium))) - 0.5)

# Metadata:
split_colnames = strsplit(colnames(tumour_medium), '_')
samples = c()
cell_types = c()
individuals = c()
for(row_split in split_colnames){
  individuals = c(individuals, row_split[1])
  samples = c(samples, row_split[2])
  cell_types = c(cell_types, row_split[3])
}
metadata = data.frame(individual=individuals, sample=samples, cell_type=cell_types,
                      row.names=colnames(tumour_medium))
more_meta = read.csv('./0Data/scRNAseq/CRC_atlas/expression_data/metadata.csv', row.names=1)
states = c()
for(samp in metadata$sample){
  states = c(states, more_meta[samp, 'Sample.Source'])
}
metadata$state = states





# --------------------------------------------------------------------------
# --- Do models' fluxes change between using a normal and tumour medium? ---
# --------------------------------------------------------------------------

# 1. Perform the statistical tests:
test_results_nvst = data.frame(p_value=rep(NA, dim(tumour_medium)[2]),
p_adjust=rep(NA, dim(tumour_medium)[2]))
row.names(test_results_nvst) = colnames(tumour_medium)
pb <- txtProgressBar(min=0, max=dim(tumour_medium)[2], style=3)
i=0
for(model in colnames(tumour_medium)){
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
write.csv(test_results_nvst, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/4_tumourMedium/normal_vs_tumour.csv')

# 2. Percentage of models that differ
sum(test_results_nvst$p_adjust < 0.05) / 172 * 100 # 71 models differ significantly, i.e., ~41%

# 3. Percentage of models that differ, per cell-type
x = table(metadata[rownames(test_results_nvst), 'cell_type'][test_results_nvst$p_adjust < 0.05]) /
  table(metadata$cell_type) * 100

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
  ggplot2::theme(axis.text.x=ggplot2::element_blank()) +
  #ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=7)) +
  ggplot2::ylab('Number of Models') + ggplot2::xlab('Cell Types')

ggplot2::ggplot(perc_df2, ggplot2::aes(x=cell_type, y=difference, fill=cell_type)) +
  ggplot2::geom_bar(stat='identity') + ggplot2::scale_fill_manual(values=ct_colors) +
  ggplot2::theme(axis.text.x=ggplot2::element_blank()) +
  ggplot2::ylab('Percentage') + ggplot2::xlab('Cell Types') + ggplot2::ylim(0, 60)

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
diff_models = colnames(tumour_medium)[test_results_nvst$p_adjust < 0.05]

# 1. For each model that differs, evaluate the differences:
reaction_differences = tumour_medium - normal_medium[, colnames(tumour_medium)]
write.csv(reaction_differences,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/4_tumourMedium/rxn_differences.csv')
reaction_FC = tumour_medium / normal_medium[, colnames(tumour_medium)]
write.csv(reaction_FC,
          './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/4_tumourMedium/rxn_FC.csv')

# 2. Metadata of only the differing reactions:
metadata_diff = metadata
metadata_diff$if_diff = rep('No', dim(metadata)[1])
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
                     './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/4_tumourMedium/rxn_changes.json')
# 3.3. Get % of reactions in subsystems that changed:
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
# 3.3.2. Mean percentage per pathway across all models:
mean_perc_models = rowMeans(perc_changed_paths)
names(mean_perc_models)[mean_perc_models == 0] # 2 pathways are not affected at all in any model, which
                                              # are dietary fiber binding, and CoA metabolism
top15_most = sort(mean_perc_models, decreasing=T)[1:15] # Top 20 most affected pathways
top15_least = sort(mean_perc_models)[1:15] # Top 20 least affected pathways
# 3.3. Visualize the percentages:
stt = metadata_diff$state
stt[stt!='Normal Matched'] = 'Tumour'
state_colors2 = c('#52854c', '#b30000')
names(state_colors2) = c('Normal Matched', 'Tumour')
heat_perc = ComplexHeatmap::Heatmap(perc_changed_paths[c(names(top15_least), names(top15_most)),
                                                       rownames(metadata_diff)],
                                  column_split = metadata_diff$cell_type,
                                  row_split = c(rep('Least', 15), rep('Most', 15)),
                                  column_title = NULL, row_title = NULL,
                                  show_column_names = FALSE,
                                  row_names_gp = grid::gpar(fontsize=7),
                                  heatmap_legend_param = list(title='% Affected reactions',
                                                              title_position = "lefttop-rot"),
                                  top_annotation = ComplexHeatmap::HeatmapAnnotation(
                                    cell_type=metadata_diff$cell_type,
                                    state=stt,
                                    if_diff=metadata_diff$if_diff,
                                    col=list(cell_type=ct_colors, state=state_colors2,
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
# 3.4. Distribution of the types of change in the top 15 subsystems:
df_plot = c()
for(model in names(rxn_changes)){
  model_df = rxn_changes[[model]]
  temp_df = model_df[model_df$state!='Unchanged' & model_df$pathway%in%names(top15_most),
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



