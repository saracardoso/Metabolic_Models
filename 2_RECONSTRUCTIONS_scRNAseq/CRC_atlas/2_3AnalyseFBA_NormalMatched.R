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

important_rxns = jsonlite::read_json('./GENERAL/utility_data/important_reactions_Tcells.json',
                                     simplifyVector=TRUE)



# Read FBA data:
fba_fluxes = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/normal_FBA.csv',
                      row.names=1, check.names = FALSE)

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

# Sigmoid Transform fluxes:
fba_fluxes = 2 * ((1 / (1 + exp(-fba_fluxes))) - 0.5)

# Filter reactions not present (zero flux) in any model
fba_fluxes_filtered = fba_fluxes[rowSums(fba_fluxes)!=0,]

# --------------
# - Clustering -
# --------------

# Perform PCA with all reactions:
pca_res = irlba::irlba(A=t(fba_fluxes_filtered), nv=50)
rxn_embeddings = pca_res$u %*% diag(pca_res$d)

pca_plot_df = cbind(rxn_embeddings[, 1:2], metadata)
colnames(pca_plot_df)[1:2] = c('PC1', 'PC2')
plt = ggplot2::ggplot(pca_plot_df, ggplot2::aes(PC1, PC2))

# Plot coloring by individuals, samples, states, or cell types
plt + ggplot2::geom_point(ggplot2::aes_string(colour='individual'), size=2) +
  ggplot2::scale_colour_manual(values=individual_colors) +
  ggplot2::labs(colour='Individual')
plt + ggplot2::geom_point(ggplot2::aes_string(colour='samples'), size=2) +
  ggplot2::theme(legend.position='none')

plt + ggplot2::geom_point(ggplot2::aes_string(colour='state'), size=2) +
  ggplot2::scale_colour_manual(values=state_colors) +
  ggplot2::labs(colour='State')

plt + ggplot2::geom_point(ggplot2::aes_string(colour='cell_type'), size=2) +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::labs(colour='Cell Type')

# Simplify cell-type annotations to see if we can see some sort of separation
metadata$simple_cts = metadata$cell_type
metadata$simple_cts[metadata$cell_type%in%c('Naive CD4 Tcells',
                                            'Naive CD8 Tcells')] = 'Naive'
metadata$simple_cts[metadata$cell_type%in%c('Proliferative CD4 Tcells',
                                            'Proliferative CD8 Tcells')] = 'Proliferative'
metadata$simple_cts[metadata$cell_type%in%c('Memory CD4 Tcells',
                                            'Memory CD8 Tcells')] = 'Memory'
metadata$simple_cts[metadata$cell_type%in%c('Follicular CD4 Tcells',
                                            'IL17+ CD8 Tcells')] = 'Effector CD4'
simple_ct_colors = ct_colors[c(1,3,5,7,9,10)]
names(simple_ct_colors)[2:5] = c('Effector CD4', 'Memory', 'Naive', 'Proliferative')
pca_plot_df = cbind(rxn_embeddings[, 1:5], metadata)
colnames(pca_plot_df)[1:5] = c('PC1', 'PC2', 'PC3', 'PC4', 'PC5')

plt = ggplot2::ggplot(pca_plot_df, ggplot2::aes(PC1, PC2, colour=simple_cts))
plt + ggplot2::geom_point(size=2) +
  ggplot2::scale_colour_manual(values=simple_ct_colors) +
  ggplot2::labs(colour='Cell Type')

plt = ggplot2::ggplot(pca_plot_df[metadata$simple_cts%in%c('Cytotoxic CD8 Tcells', 'Memory', 'Regulatory CD4 Tcells'),],
                      ggplot2::aes(PC1, PC2, colour=simple_cts))
plt + ggplot2::geom_point(size=2) +
  ggplot2::scale_colour_manual(values=simple_ct_colors) +
  ggplot2::labs(colour='Cell Type')





# ------------------------------------------------------
# - Differential Flux Analysis and Enrichment Analysis -
# ------------------------------------------------------

HumanGEM_subsystems_rxns = jsonlite::read_json('./GENERAL/utility_data/subsystems_reactions_mapping.json',
                                               simplifyVector = TRUE)

# 1. Reaction Set Enrichment Analysis
# For each ct, calculate fold change between ct and others. Sort fc
# decreasingly and perform EA.
rsea = list()
for(ct in unique(metadata$cell_type)){
  message(ct)
  ct_rxns_means = rowMeans(fba_fluxes_filtered[,metadata$cell_type==ct])
  other_rxns_means = rowMeans(fba_fluxes_filtered[,metadata$cell_type!=ct])
  fc = ct_rxns_means - other_rxns_means
  names(fc) = row.names(fba_fluxes_filtered)
  fc = sort(fc, decreasing = TRUE)
  rsea[[ct]] = fgsea::fgsea(HumanGEM_subsystems_rxns, fc, eps=0,
                            scoreType='pos')
}
jsonlite::write_json(rsea, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/rsea.json')
# Make dotplot (color -> enrichment score (ES); size -> -log(p.value)) -> only
# show pathways that are enriched in at least one of the cell-types
# What are the enriched paths for each ct?
interesting_pathways = c()
for(ct in names(rsea)){
  interesting_pathways = unique(c(interesting_pathways, rsea[[ct]]$pathway[rsea[[ct]]$padj<0.05 &
                                                                             rsea[[ct]]$ES>0.5]))
}
dotplot_df = c()
for(ct in names(rsea)){
  dotplot_df = rbind(dotplot_df,
                     cbind(rsea[[ct]][rsea[[ct]]$pathway%in%interesting_pathways,
                                      c('pathway', 'padj', 'ES')],
                           rep(ct, length(interesting_pathways))))
}
colnames(dotplot_df)[4] = 'cell_type'
dotplot_df$pathway = factor(dotplot_df$pathway,
                            levels=sort(unique(dotplot_df$pathway), decreasing=TRUE))
dotplot_df$padj = -log10(dotplot_df$padj)
dotplot_df$padj_discrete = dotplot_df$padj
dotplot_df$padj_discrete[dotplot_df$padj_discrete <= (-log10(0.05))] = 0
dotplot_df$padj_discrete[dotplot_df$padj_discrete > (-log10(0.05))] = 1
ggplot2::ggplot(dotplot_df, ggplot2::aes(x=cell_type, y=pathway, colour=ES,
                                         size=padj_discrete)) +
  ggplot2::geom_point() +
  ggplot2::scale_colour_gradient2(low='blue', mid='white', high='red',
                                  midpoint=0.5) +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1)) +
  ggplot2::scale_size(breaks=c(0, 1), labels=c("> 0.05","< 0.05"), guide="legend") +
  ggplot2::labs(colour='Enrichment Score', size='Adjusted p-value') +
  ggplot2::xlab('') + ggplot2::ylab('')


# 2. Differential reaction flux
# What are the differential fluxes for each ct? (only consider those with fc
# positive)
diff_expressed_rxns = c()
for(ct in unique(metadata$cell_type)){
  message('-',ct)
  # - Only consider reactions whose fold change is positive
  ct_rxns_means = rowMeans(fba_fluxes_filtered[,metadata$cell_type==ct])
  other_rxns_means = rowMeans(fba_fluxes_filtered[,metadata$cell_type!=ct])
  fc = ct_rxns_means - other_rxns_means
  names(fc) = row.names(fba_fluxes_filtered)
  rxns_toTest = names(fc)[fc>0]
  # - For those reactions, find the ones that are significantly different
  test_meta1 = metadata$cell_type
  test_meta1 = data.frame(cell_type=test_meta1,
                          row.names=colnames(fba_fluxes_filtered))
  test_meta1$cell_type[test_meta1$cell_type!=ct] = 'other'
  p_val = sapply(rxns_toTest,
                 FUN = function(x) {
                   return(wilcox.test(as.numeric(fba_fluxes_filtered[x, ]) ~ as.vector(test_meta1[, "cell_type"]))$p.value)
                 }
  )
  pval_adjust = p.adjust(p_val, method='fdr')
  res_ct = data.frame(reaction=rxns_toTest, pval_adjust,
                      pval=p_val, cell_type=ct,
                      row.names=paste(rxns_toTest, ct, sep='_'))
  diff_expressed_rxns = rbind(diff_expressed_rxns, res_ct)
}
diff_expressed_rxns_final = diff_expressed_rxns[diff_expressed_rxns$pval_adjust<0.05,]
different_rxns = unique(diff_expressed_rxns_final$reaction)
write.csv(diff_expressed_rxns_final, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/dre.csv')
# Perform Over-representation analysis -> similar to previous enrichment?
paths_rxns_ora = c()
for(path in names(HumanGEM_subsystems_rxns)){
  path_rxns = cbind(rep(path, length(HumanGEM_subsystems_rxns[[path]])),
                    HumanGEM_subsystems_rxns[[path]])
  if(length(paths_rxns_ora)==0) paths_rxns_ora = path_rxns
  else paths_rxns_ora = rbind(paths_rxns_ora, path_rxns)
}
ora_res = list()
for(ct in unique(metadata$cell_type)){
  message(ct)
  ct_diff_rxns = diff_expressed_rxns_final$reaction[diff_expressed_rxns_final$pval_adjust<0.05 &
                                                      diff_expressed_rxns_final$cell_type==ct]
  temp = clusterProfiler::enricher(gene=ct_diff_rxns,
                                   universe=rownames(fba_fluxes_filtered),
                                   TERM2GENE=paths_rxns_ora)
  if(!is.null(temp)) ora_res[[ct]] = temp@result
}
jsonlite::write_json(ora_res, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/1_control_analysis/FBA/ora.json')
# Make dotplot (color -> Reaction ratio; size -> -log(p.value)) -> only
# show pathways that are enriched in at least one of the cell-types
# What are the enriched paths for each ct?
interesting_pathways = c()
for(ct in names(ora_res)){
  interesting_pathways = unique(c(interesting_pathways, rownames(ora_res[[ct]])[ora_res[[ct]]$p.adjust<0.05]))
}
dotplot_df = c()
for(ct in names(ora_res)){
  path_info = ora_res[[ct]][rownames(ora_res[[ct]])%in%interesting_pathways,
                            c('ID', 'GeneRatio', 'p.adjust')]
  dotplot_df = rbind(dotplot_df, cbind(path_info, rep(ct, dim(path_info)[1])))
}
dotplot_df$GeneRatio = gsubfn::strapply(dotplot_df$GeneRatio,
                                        "([0-9]+)/([0-9]+)",
                                        ~ as.numeric(x) / as.numeric(y),
                                        backref = -2, simplify = c)
colnames(dotplot_df)[2] = 'ReactionRatio'
colnames(dotplot_df)[4] = 'cell_type'
dotplot_df$ID = factor(dotplot_df$ID,
                       levels=sort(unique(dotplot_df$ID), decreasing=TRUE))
dotplot_df$p.adjust = -log10(dotplot_df$p.adjust)
dotplot_df$padj_discrete = dotplot_df$p.adjust
dotplot_df$padj_discrete[dotplot_df$padj_discrete <= (-log10(0.05))] = 0
dotplot_df$padj_discrete[dotplot_df$padj_discrete > (-log10(0.05))] = 1
ggplot2::ggplot(dotplot_df, ggplot2::aes(x=cell_type, y=ID,
                                         colour=ReactionRatio,
                                         size=padj_discrete)) +
  ggplot2::geom_point() +
  ggplot2::scale_colour_gradient2(low='blue', mid='white', high='red',
                                  midpoint=0.5) +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1)) +
  ggplot2::scale_size(breaks=c(0, 1), labels=c("> 0.05","< 0.05"), guide="legend") +
  ggplot2::labs(colour='Reaction Ratio', size='Adjusted p-value') +
  ggplot2::xlab('') + ggplot2::ylab('')


# 3. New PCA with only differential reactions -> does it improve anything?
pca_res = irlba::irlba(A=t(fba_fluxes_filtered[different_rxns,]), nv=50)
rxn_embeddings = pca_res$u %*% diag(pca_res$d)

pca_plot_df = cbind(rxn_embeddings[, 1:2], metadata)
colnames(pca_plot_df)[1:2] = c('PC1', 'PC2')
plt = ggplot2::ggplot(pca_plot_df, ggplot2::aes(PC1, PC2))
# Plot coloring by individuals, samples, states, or cell types
plt + ggplot2::geom_point(ggplot2::aes_string(colour='individual'), size=2) +
  ggplot2::scale_colour_manual(values=individual_colors) +
  ggplot2::labs(colour='Individual')
plt + ggplot2::geom_point(ggplot2::aes_string(colour='samples'), size=2) +
  ggplot2::theme(legend.position='none')

plt + ggplot2::geom_point(ggplot2::aes_string(colour='state'), size=2) +
  ggplot2::scale_colour_manual(values=state_colors) +
  ggplot2::labs(colour='State')

plt + ggplot2::geom_point(ggplot2::aes_string(colour='cell_type'), size=2) +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::labs(colour='Cell Type')
# Simplify cell-type annotations to see if we can see some sort of separation
metadata$simple_cts = metadata$cell_type
metadata$simple_cts[metadata$cell_type%in%c('Naive CD4 Tcells',
                                            'Naive CD8 Tcells')] = 'Naive'
metadata$simple_cts[metadata$cell_type%in%c('Proliferative CD4 Tcells',
                                            'Proliferative CD8 Tcells')] = 'Proliferative'
metadata$simple_cts[metadata$cell_type%in%c('Memory CD4 Tcells',
                                            'Memory CD8 Tcells')] = 'Memory'
metadata$simple_cts[metadata$cell_type%in%c('Follicular CD4 Tcells',
                                            'IL17+ CD8 Tcells')] = 'Effector CD4'
simple_ct_colors = ct_colors[c(1,3,5,7,9,10)]
names(simple_ct_colors)[2:5] = c('Effector CD4', 'Memory', 'Naive', 'Proliferative')
pca_plot_df = cbind(rxn_embeddings[, 1:5], metadata)
colnames(pca_plot_df)[1:5] = c('PC1', 'PC2', 'PC3', 'PC4', 'PC5')

plt = ggplot2::ggplot(pca_plot_df, ggplot2::aes(PC1, PC2, colour=simple_cts))
plt + ggplot2::geom_point(size=2) +
  ggplot2::scale_colour_manual(values=simple_ct_colors) +
  ggplot2::labs(colour='Cell Type')
