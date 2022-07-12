ct_colors = c('#b30000', '#999999', '#cc79a7', '#b3c9e6', '#4e84c4', '#c3d7a4', '#52854c', '#fdb981',
              '#d16103', '#293352')
names(ct_colors) = c('Cytotoxic CD8 Tcells', 'Follicular CD4 Tcells', 'IL17+ CD4 Tcells',
                     'Memory CD4 Tcells', 'Memory CD8 Tcells', 'Naive CD4 Tcells', 'Naive CD8 Tcells',
                     'Proliferative CD4 Tcells', 'Proliferative CD8 Tcells',  'Regulatory CD4 Tcells')
state_colors = c('#52854c', '#b30000')
names(state_colors) = c('Normal Matched', 'Tumour')

code_dir = './GENERAL/code/R'
for(file in list.files(code_dir, full.names=TRUE)) source(file)





# -----------------
# --- Load Data ---
# -----------------

metadata = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/metadata.csv', row.names=1)
fba_normalBlood_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/Normal_Blood.csv',
                               row.names=1, check.names = FALSE)
fba_normalBlood = 2 * ((1 / (1 + exp(-fba_normalBlood_raw))) - 0.5)




# ---------------
# --- BIOMASS ---
# ---------------

plot_df = data.frame(flux=as.numeric(fba_normalBlood['MAR13082', rownames(metadata)]),
                     ct=metadata$cell_type, state=metadata$state)
ggplot2::ggplot(plot_df, ggplot2::aes(x=ct, y=flux, colour=ct)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') + ggplot2::ylab('Flux (mmol/gDW/h)')

ggplot2::ggplot(plot_df, ggplot2::aes(x=state, y=flux, colour=ct)) +
  ggplot2::geom_boxplot(outlier.shape=NA) + ggplot2::facet_wrap(ggplot2::vars(ct)) +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') + ggplot2::ylab('Flux (mmol/gDW/h)')





# ---------------------------------
# --- Sources of FADH2 and NADH ---
# ---------------------------------

xx = source_nadh_fadh2_ct(fba_normalBlood[, rownames(metadata)], metadata)

values = c()
for(k in 1:dim(xx)[2]) values = c(values, xx[,k])
source_energy = c(rep('FADH2', dim(xx)[1]), rep('NADH', dim(xx)[1]), rep('FADH2', dim(xx)[1]),
                  rep('NADH', dim(xx)[1]*3))
pathways = c(rep('FAO', dim(xx)[1]*2), rep('TCA', dim(xx)[1]*2),
             rep('Glycolysis', dim(xx)[1]), rep('Glutaminolysis', dim(xx)[1]))
cell_types = rep(metadata[rownames(xx),'cell_type'], dim(xx)[2])

plot_df = data.frame(values, source_energy, pathways, cell_types)
fadh2_plot = plot_df[plot_df$source_energy=='FADH2',]
nadh_plot = plot_df[plot_df$source_energy=='NADH',]


ggplot2::ggplot(fadh2_plot, ggplot2::aes(x=pathways, y=values, colour=pathways)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1),
                 legend.position='none') + ggplot2::ylab('Flux')

ggplot2::ggplot(nadh_plot, ggplot2::aes(x=pathways, y=values, colour=pathways)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1),
                 legend.position='none') + ggplot2::ylab('Flux')





# ------------------
# --- FAS Uptake ---
# ------------------

fau_rxns = jsonlite::read_json('./GENERAL/utility_data/important_reactions_Tcells.json',
                                     simplifyVector=TRUE)$uptakes$`fatty acids`

fau = colSums(na.omit(fba_normalBlood[fau_rxns, rownames(metadata)]))
plot_df = data.frame(flux=fau, ct=metadata$cell_type, state=metadata$state)

ggplot2::ggplot(plot_df, ggplot2::aes(x=ct, y=flux, colour=ct)) +
  ggplot2::geom_boxplot(outlier.shape=NA) + ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1),
                 legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) + ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')

ggplot2::ggplot(plot_df, ggplot2::aes(x=state, y=flux, colour=ct)) +
  ggplot2::geom_boxplot(outlier.shape=NA) + ggplot2::geom_jitter(width = 0.15) +
  ggplot2::facet_wrap(ggplot2::vars(ct)) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1),
                 legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) + ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')





# ----------------------
# --- Medium Changes ---
# ----------------------


# 1. No tryptophan:
B2_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/medium_changes/B2.csv',
                         row.names=1, check.names=FALSE)
B2_fluxes = 2 * ((1 / (1 + exp(-B2_fluxes_raw))) - 0.5)
plot_df = data.frame(values=c(as.numeric(fba_normalBlood['MAR13082', rownames(metadata)]),
                              as.numeric(B2_fluxes['MAR13082', rownames(metadata)])),
                     cell_types=rep(metadata[,'cell_type'], 2),
                     medium=factor(c(rep('Normal', dim(metadata)[1]),
                                     rep('No Tryptophan', dim(metadata)[1])),
                                   levels=c('Normal', 'No Tryptophan')),
                     paired=c(1:196, 1:196))

ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2,
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')


# 2. No Oxygen:
B4_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/medium_changes/B4.csv',
                         row.names=1, check.names=FALSE)
B4_fluxes = 2 * ((1 / (1 + exp(-B4_fluxes_raw))) - 0.5)
plot_df = data.frame(values=c(as.numeric(fba_normalBlood['MAR13082', rownames(metadata)]),
                              as.numeric(B4_fluxes['MAR13082', rownames(metadata)])),
                     cell_types=rep(metadata[,'cell_type'], 2),
                     state=rep(metadata[,'state'], 2),
                     medium=factor(c(rep('Normal', dim(metadata)[1]),
                                     rep('No Oxygen', dim(metadata)[1])),
                                   levels=c('Normal', 'No Oxygen')),
                     paired=c(1:196, 1:196))
ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2,
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')
ggplot2::ggplot(plot_df, ggplot2::aes(x=state, y=values, colour=medium)) +
  ggplot2::geom_boxplot(outlier.shape=NA) + ggplot2::facet_wrap(ggplot2::vars(cell_types)) + 
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1), legend.position='right') +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')
# 2.1. Reaction catalised by oxygen
plot_df = data.frame(values=c(as.numeric(fba_normalBlood['MAR04776', rownames(metadata)]),
                              as.numeric(B4_fluxes['MAR04776', rownames(metadata)])),
                     cell_types=rep(metadata[,'cell_type'], 2),
                     medium=factor(c(rep('Normal', dim(metadata)[1]),
                                     rep('No Oxygen', dim(metadata)[1])),
                                   levels=c('Normal', 'No Oxygen')),
                     paired=c(1:196, 1:196))
ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2,
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')
# 2.2. Biomass difference vs oxygen uptake:
bmss_diff = abs(as.numeric(B4_fluxes['MAR13082', rownames(metadata)]) -
                  as.numeric(fba_normalBlood['MAR13082', rownames(metadata)]))
plt_df = data.frame(biomass_diff=bmss_diff,
                    uptake=(-as.numeric(fba_normalBlood['MAR09048', rownames(metadata)])))
ggplot2::ggplot(plt_df, ggplot2::aes(x=uptake, y=biomass_diff)) +
  ggplot2::geom_point() + ggplot2::geom_smooth(method='lm', se=FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::xlab('Oxygen Uptake (mmol/gDW/h)') +
  ggplot2::ylab('Absolute Biomass Difference (mmol/gDW/h)')


# 3. No Glutamine:
B82_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/medium_changes/B8.2_and_B9.csv',
                         row.names=1, check.names=FALSE)
B82_fluxes = 2 * ((1 / (1 + exp(-B82_fluxes_raw))) - 0.5)
plot_df = data.frame(values=c(as.numeric(fba_normalBlood['MAR13082', rownames(metadata)]),
                              as.numeric(B82_fluxes['MAR13082', rownames(metadata)])),
                     cell_types=rep(metadata[,'cell_type'], 2),
                     state=rep(metadata[,'state'], 2),
                     medium=factor(c(rep('Normal', dim(metadata)[1]),
                                     rep('No Glutamine', dim(metadata)[1])),
                                   levels=c('Normal', 'No Glutamine')),
                     paired=c(1:196, 1:196))
ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2,
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')
# 3.1. DNA production
plot_df = data.frame(values=c(as.numeric(fba_normalBlood['MAR07160', rownames(metadata)]),
                              as.numeric(B82_fluxes['MAR07160', rownames(metadata)])),
                     cell_types=rep(metadata[,'cell_type'], 2),
                     state=rep(metadata[,'state'], 2),
                     medium=factor(c(rep('Normal', dim(metadata)[1]),
                                     rep('No Glutamine', dim(metadata)[1])),
                                   levels=c('Normal', 'No Glutamine')),
                     paired=c(1:196, 1:196))
ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2,
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')
# 3.2. RNA production
plot_df = data.frame(values=c(as.numeric(colSums(fba_normalBlood[c('MAR07161', 'MAR07162'),
                                                            colnames(B82_fluxes)])),
                              as.numeric(colSums(B82_fluxes[c('MAR07161', 'MAR07162'),
                                                            colnames(B82_fluxes)]))),
                     cell_types=rep(metadata[,'cell_type'], 2),
                     state=rep(metadata[,'state'], 2),
                     medium=factor(c(rep('Normal', dim(metadata)[1]),
                                     rep('No Glutamine', dim(metadata)[1])),
                                   levels=c('Normal', 'No Glutamine')),
                     paired=c(1:196, 1:196))
ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2,
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')


# 4. No Nucleotides:
# 4.1. Biomass
B17_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/medium_changes/B17.csv',
                          row.names=1, check.names=FALSE)
B17_fluxes = 2 * ((1 / (1 + exp(-B17_fluxes_raw))) - 0.5)
plot_df = data.frame(values=c(as.numeric(fba_normalBlood['MAR13082', rownames(metadata)]),
                              as.numeric(B17_fluxes['MAR13082', rownames(metadata)])),
                     cell_types=rep(metadata[,'cell_type'], 2),
                     state=rep(metadata[,'state'], 2),
                     medium=factor(c(rep('Normal', dim(metadata)[1]),
                                     rep('No Nucleotides', dim(metadata)[1])),
                                   levels=c('Normal', 'No Nucleotides')),
                     paired=c(1:196, 1:196))
ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2,
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')
# 4.2. Production of nucleotides
plot_df = data.frame(values=c(as.numeric(fba_normalBlood['MAR07160', rownames(metadata)]),
                              as.numeric(B17_fluxes['MAR07160', rownames(metadata)])),
                     cell_types=rep(metadata[,'cell_type'], 2),
                     state=rep(metadata[,'state'], 2),
                     medium=factor(c(rep('Normal', dim(metadata)[1]),
                                     rep('No Nucleotides', dim(metadata)[1])),
                                   levels=c('Normal', 'No Nucleotides')),
                     paired=c(1:196, 1:196))
ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2,
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')


