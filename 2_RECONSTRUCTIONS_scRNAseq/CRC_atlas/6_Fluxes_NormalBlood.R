ct_colors = c('#b30000', '#999999', '#cc79a7', '#b3c9e6', '#4e84c4', '#c3d7a4', '#52854c', '#fdb981',
              '#d16103', '#293352')
names(ct_colors) = c('Cytotoxic CD8 Tcells', 'Follicular CD4 Tcells', 'IL17+ CD4 Tcells',
                     'Memory CD4 Tcells', 'Memory CD8 Tcells', 'Naive CD4 Tcells', 'Naive CD8 Tcells',
                     'Proliferative CD4 Tcells', 'Proliferative CD8 Tcells',  'Regulatory CD4 Tcells')
state_colors = c('#52854c', '#b30000')
names(state_colors) = c('Normal Matched', 'Tumour')
cms_colors = c('#d4cfcf', '#cc79a7', '#b3c9e6', '#c3d7a4', '#52854c', '#fdb981')
names(cms_colors) = c('CMS3', 'Mixed', 'CMS1', 'CMS2', 'Normal Matched', 'CMS4')

#code_dir = './GENERAL/code/R'
#for(file in list.files(code_dir, full.names=TRUE)) source(file)
source('./GENERAL/code/R/bar_graphs.R')





# -----------------
# --- Load Data ---
# -----------------

metadata = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/metadata.csv', row.names=1)
pfba_normalBlood_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/Normal_Blood_pFBA.csv',
                               row.names=1, check.names = FALSE)
pfba_normalBlood = 2 * ((1 / (1 + exp(-pfba_normalBlood_raw))) - 0.5)

pfba_normalBlood_biomass_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/Normal_Blood_pFBA_biomass.csv',
                               row.names=1, check.names = FALSE)
pfba_normalBlood_biomass = 2 * ((1 / (1 + exp(-pfba_normalBlood_biomass_raw))) - 0.5)




# ---------------
# --- BIOMASS ---
# ---------------

plot_df = data.frame(flux=as.numeric(pfba_normalBlood['MAR13082', ]),
                     ct=metadata[colnames(pfba_normalBlood), 'cell_type'],
                     state=metadata[colnames(pfba_normalBlood), 'state'],
                     CMS=metadata[colnames(pfba_normalBlood), 'CMS'])
# OR:
plot_df = data.frame(flux=as.numeric(pfba_normalBlood_biomass['MAR13082', ]),
                     ct=metadata[colnames(pfba_normalBlood_biomass), 'cell_type'],
                     state=metadata[colnames(pfba_normalBlood_biomass), 'state'],
                     CMS=metadata[colnames(pfba_normalBlood_biomass), 'CMS'])

ggplot2::ggplot(plot_df, ggplot2::aes(x=ct, y=flux, colour=ct)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') + ggplot2::ylab('Flux (mmol/gDW/h)')

ggplot2::ggplot(plot_df, ggplot2::aes(x=state, y=flux, colour=ct)) +
  ggplot2::geom_boxplot(outlier.shape=NA) + ggplot2::facet_wrap(ggplot2::vars(ct)) +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=7), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') + ggplot2::ylab('Flux (mmol/gDW/h)')

ggplot2::ggplot(plot_df, ggplot2::aes(x=CMS, y=flux, colour=ct)) +
  ggplot2::geom_boxplot(outlier.shape=NA) + ggplot2::facet_wrap(ggplot2::vars(ct)) +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=7), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') + ggplot2::ylab('Flux (mmol/gDW/h)')





# ----------------------
# --- ATP production ---
# ----------------------

plot_df = data.frame(flux=as.numeric(pfba_normalBlood['MAR06916', ]),
                     ct=metadata[colnames(pfba_normalBlood), 'cell_type'],
                     state=metadata[colnames(pfba_normalBlood), 'state'],
                     CMS=metadata[colnames(pfba_normalBlood), 'CMS'])
# OR:
plot_df = data.frame(flux=as.numeric(pfba_normalBlood_biomass['MAR06916', ]),
                     ct=metadata[colnames(pfba_normalBlood_biomass), 'cell_type'],
                     state=metadata[colnames(pfba_normalBlood_biomass), 'state'],
                     CMS=metadata[colnames(pfba_normalBlood_biomass), 'CMS'])

ggplot2::ggplot(plot_df, ggplot2::aes(x=ct, y=flux, colour=ct)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') + ggplot2::ylab('Flux (mmol/gDW/h)')

ggplot2::ggplot(plot_df, ggplot2::aes(x=state, y=flux, colour=ct)) +
  ggplot2::geom_boxplot(outlier.shape=NA) + ggplot2::facet_wrap(ggplot2::vars(ct)) +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=7), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') + ggplot2::ylab('Flux (mmol/gDW/h)')

ggplot2::ggplot(plot_df, ggplot2::aes(x=CMS, y=flux, colour=ct)) +
  ggplot2::geom_boxplot(outlier.shape=NA) + ggplot2::facet_wrap(ggplot2::vars(ct)) +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=7), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') + ggplot2::ylab('Flux (mmol/gDW/h)')





# ---------------------------------
# --- Biomass vs ATP production ---
# ---------------------------------

plot_df = data.frame(atp=as.numeric(pfba_normalBlood['MAR06916', ]),
                     biomass=as.numeric(pfba_normalBlood['MAR13082', ]),
                     ct=metadata[colnames(pfba_normalBlood), 'cell_type'],
                     state=metadata[colnames(pfba_normalBlood), 'state'],
                     CMS=metadata[colnames(pfba_normalBlood), 'CMS'])

ggplot2::ggplot(plot_df, ggplot2::aes(x=atp, y=biomass, colour=ct)) +
  ggplot2::geom_point(size=2) +
  ggplot2::geom_vline(xintercept=0.5, color="#404040", linetype='dashed') +
  ggplot2::geom_hline(yintercept=0.005, color="#404040", linetype='dashed') +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('ATP production Flux (mmol/gDW/h)') + ggplot2::ylab('Biomass Flux (mmol/gDW/h)')

# Proliferative models with relatively high biomass and low ATP: 13 out of 19.
sum(plot_df$atp < 0.5 & plot_df$biomass > 0.005 & plot_df$ct%in%c('Proliferative CD4 Tcells', 'Proliferative CD8 Tcells'))
sum(plot_df$ct%in%c('Proliferative CD4 Tcells', 'Proliferative CD8 Tcells'))

sum(!plot_df$ct%in%c('Proliferative CD4 Tcells', 'Proliferative CD8 Tcells')) # 123
# Non-proliferative models with relatively high biomass and low ATP: 28
sum(plot_df$atp < 0.5 & plot_df$biomass > 0.005 & !plot_df$ct%in%c('Proliferative CD4 Tcells', 'Proliferative CD8 Tcells'))
# Non-proliferative models with relatively low biomass and high ATP: 73
sum(plot_df$atp > 0.5 & plot_df$biomass < 0.005 & !plot_df$ct%in%c('Proliferative CD4 Tcells', 'Proliferative CD8 Tcells'))
# Non-proliferative models with relatively low biomass and low ATP: 7
sum(plot_df$atp < 0.5 & plot_df$biomass < 0.005 & !plot_df$ct%in%c('Proliferative CD4 Tcells', 'Proliferative CD8 Tcells'))
# Non-proliferative models with relatively high biomass and high ATP: 15
sum(plot_df$atp > 0.5 & plot_df$biomass > 0.005 & !plot_df$ct%in%c('Proliferative CD4 Tcells', 'Proliferative CD8 Tcells'))

ggplot2::ggplot(plot_df, ggplot2::aes(x=atp, y=biomass, colour=state)) +
  ggplot2::geom_point(size=2) +
  ggplot2::geom_vline(xintercept=0.5, color="#404040", linetype='dashed') +
  ggplot2::geom_hline(yintercept=0.005, color="#404040", linetype='dashed') +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::scale_colour_manual(values=state_colors) +
  ggplot2::xlab('ATP production Flux (mmol/gDW/h)') + ggplot2::ylab('Biomass Flux (mmol/gDW/h)')

# Models with high atp and high biomass:
sum(plot_df$atp > 0.5 & plot_df$biomass > 0.005) # 15 (all are non-proliferative cell-types)
# Tumour vs normal models for those with high atp and biomass:
table(plot_df$state[plot_df$atp > 0.5 & plot_df$biomass > 0.005]) # 10 are tumour and 5 are normal

# Models with low atp and low biomass:
sum(plot_df$atp < 0.5 & plot_df$biomass < 0.005) # 12 (7 are non-proliferative cell-types)
# Tumour vs normal models for those with low atp and biomass:
table(plot_df$state[plot_df$atp < 0.5 & plot_df$biomass < 0.005]) # 9 are tumour and 3 are normal





# ---------------------------------
# --- Sources of FADH2 and NADH ---
# ---------------------------------

xx = source_nadh_fadh2_ct(pfba_normalBlood, metadata[colnames(pfba_normalBlood), ])
# OR:
xx = source_nadh_fadh2_ct(pfba_normalBlood_biomass, metadata[colnames(pfba_normalBlood_biomass), ])

values = c()
for(k in 1:dim(xx)[2]) values = c(values, xx[,k])
source_energy = c(rep('FADH2', dim(xx)[1]), rep('NADH', dim(xx)[1]), rep('FADH2', dim(xx)[1]),
                  rep('NADH', dim(xx)[1]*3))
pathways = c(rep('FAO', dim(xx)[1]*2), rep('TCA', dim(xx)[1]*2),
             rep('Glycolysis', dim(xx)[1]), rep('Glutaminolysis', dim(xx)[1]))
cell_types = rep(metadata[rownames(xx),'cell_type'], dim(xx)[2])
samp_names = rep(rownames(xx), dim(xx)[2])

plot_df = data.frame(values, source_energy, pathways, cell_types, samp_names)
fadh2_plot = plot_df[plot_df$source_energy=='FADH2',]
nadh_plot = plot_df[plot_df$source_energy=='NADH',]

# FADH2:
ggplot2::ggplot(fadh2_plot, ggplot2::aes(x=pathways, y=values, colour=pathways)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types))+#, scales='free') +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8),
                 legend.position='none') + ggplot2::ylab('Flux')
# Only with proliferative T-cells:
ggplot2::ggplot(fadh2_plot[fadh2_plot$cell_types%in%c('Proliferative CD4 Tcells', 'Proliferative CD8 Tcells'),],
                ggplot2::aes(x=pathways, y=values, colour=pathways)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types))+#, scales='free') +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8),
                 legend.position='none') + ggplot2::ylab('Flux')

values_added_up = c()
samps_added_up = c()
cts_added_up = c()
for(ct in unique(fadh2_plot$cell_types)){
  samps_ct = unique(fadh2_plot$samp_names[fadh2_plot$cell_types==ct])
  for(samp in samps_ct){
    values_added_up = c(values_added_up, sum(fadh2_plot$values[fadh2_plot$cell_types==ct & fadh2_plot$samp_names==samp]))
    samps_added_up = c(samps_added_up, samp)
    cts_added_up = c(cts_added_up, ct)
  }
}
fadh2_plot_added_up = data.frame(values=values_added_up, cell_types=cts_added_up, samps=samps_added_up, row.names=samps_added_up)
ggplot2::ggplot(fadh2_plot_added_up, ggplot2::aes(x=cell_types, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8),
                 legend.position='none') + ggplot2::ylab('Flux')

# NADH:
ggplot2::ggplot(nadh_plot, ggplot2::aes(x=pathways, y=values, colour=pathways)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types))+#, scales='free') +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8),
                 legend.position='none') + ggplot2::ylab('Flux')
# Only with proliferative T-cells:
ggplot2::ggplot(nadh_plot[nadh_plot$cell_types%in%c('Proliferative CD4 Tcells', 'Proliferative CD8 Tcells'),],
                ggplot2::aes(x=pathways, y=values, colour=pathways)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types))+#, scales='free') +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8),
                 legend.position='none') + ggplot2::ylab('Flux')

ggplot2::ggplot(nadh_plot, ggplot2::aes(x=cell_types, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::facet_wrap(ggplot2::vars(pathways), scales='free') +
  ggplot2::scale_color_manual(values=ct_colors) +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=6),
                 legend.position='none') + ggplot2::ylab('Flux')

values_added_up = c()
samps_added_up = c()
cts_added_up = c()
for(ct in unique(nadh_plot$cell_types)){
  samps_ct = unique(nadh_plot$samp_names[nadh_plot$cell_types==ct])
  for(samp in samps_ct){
    values_added_up = c(values_added_up, sum(nadh_plot$values[nadh_plot$cell_types==ct & nadh_plot$samp_names==samp]))
    samps_added_up = c(samps_added_up, samp)
    cts_added_up = c(cts_added_up, ct)
  }
}
nadh_plot_added_up = data.frame(values=values_added_up, cell_types=cts_added_up, samps=samps_added_up, row.names=samps_added_up)
ggplot2::ggplot(nadh_plot_added_up, ggplot2::aes(x=cell_types, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_jitter(width = 0.15) +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8),
                 legend.position='none') + ggplot2::ylab('Flux')





# ------------------
# --- FAS Uptake ---
# ------------------

fau_rxns = jsonlite::read_json('./GENERAL/utility_data/important_reactions_Tcells.json',
                                     simplifyVector=TRUE)$uptakes$`fatty acids`

fau = colSums(na.omit(pfba_normalBlood[fau_rxns, ]))
plot_df = data.frame(flux=fau, ct=metadata[colnames(pfba_normalBlood), 'cell_type'],
                     state=metadata[colnames(pfba_normalBlood), 'state'],
                     cms=metadata[colnames(pfba_normalBlood), 'CMS'])

ggplot2::ggplot(plot_df, ggplot2::aes(x=ct, y=flux, colour=ct)) +
  ggplot2::geom_boxplot(outlier.shape=NA) + ggplot2::geom_jitter(width = 0.15) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8),
                 legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) + ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')

ggplot2::ggplot(plot_df, ggplot2::aes(x=state, y=flux, colour=ct)) +
  ggplot2::geom_boxplot(outlier.shape=NA) + ggplot2::geom_jitter(width = 0.15) +
  ggplot2::facet_wrap(ggplot2::vars(ct)) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8),
                 legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) + ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')

ggplot2::ggplot(plot_df, ggplot2::aes(x=cms, y=flux, colour=ct)) +
  ggplot2::geom_boxplot(outlier.shape=NA) + ggplot2::geom_jitter(width = 0.15) +
  ggplot2::facet_wrap(ggplot2::vars(ct)) +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55, hjust=1, vjust=1, size=8),
                 legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) + ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')





# ----------------------
# --- Medium Changes ---
# ----------------------


# 1. No tryptophan:
B2_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/medium_changes/B2_biomass.csv',
                         row.names=1, check.names=FALSE)
B2_fluxes = 2 * ((1 / (1 + exp(-B2_fluxes_raw))) - 0.5)
plot_df = data.frame(values=c(as.numeric(pfba_normalBlood_biomass['MAR13082', colnames(B2_fluxes)]),
                              as.numeric(B2_fluxes['MAR13082', colnames(B2_fluxes)])),
                     cell_types=rep(metadata[colnames(B2_fluxes),'cell_type'], 2),
                     medium=factor(c(rep('Normal', dim(B2_fluxes)[2]),
                                     rep('No Tryptophan', dim(B2_fluxes)[2])),
                                   levels=c('Normal', 'No Tryptophan')),
                     paired=c(1:174, 1:174))

ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey',
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')


# 2. No Oxygen:
B4_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/medium_changes/B4_biomass.csv',
                         row.names=1, check.names=FALSE)
B4_fluxes = 2 * ((1 / (1 + exp(-B4_fluxes_raw))) - 0.5)
plot_df = data.frame(values=c(as.numeric(pfba_normalBlood_biomass['MAR13082', colnames(B4_fluxes)]),
                              as.numeric(B4_fluxes['MAR13082', colnames(B4_fluxes)])),
                     cell_types=rep(metadata[colnames(B4_fluxes),'cell_type'], 2),
                     state=rep(metadata[colnames(B4_fluxes), 'state'], 2),
                     medium=factor(c(rep('Normal', dim(B4_fluxes)[2]),
                                     rep('No Oxygen', dim(B4_fluxes)[2])),
                                   levels=c('Normal', 'No Oxygen')),
                     paired=c(1:174, 1:174))

ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey',
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
plot_df = data.frame(values=c(as.numeric(pfba_normalBlood_biomass['MAR04776', colnames(B4_fluxes)]),
                              as.numeric(B4_fluxes['MAR04776', colnames(B4_fluxes)])),
                     cell_types=rep(metadata[colnames(B4_fluxes),'cell_type'], 2),
                     medium=factor(c(rep('Normal', dim(B4_fluxes)[2]),
                                     rep('No Oxygen', dim(B4_fluxes)[2])),
                                   levels=c('Normal', 'No Oxygen')),
                     paired=c(1:174, 1:174))
ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey',
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')

# 2.2. Biomass difference vs oxygen uptake:
bmss_diff = abs(as.numeric(B4_fluxes['MAR13082', colnames(B4_fluxes)]) -
                  as.numeric(pfba_normalBlood_biomass['MAR13082', colnames(B4_fluxes)]))
plt_df = data.frame(biomass_diff=bmss_diff,
                    uptake=(-as.numeric(pfba_normalBlood_biomass['MAR09048', colnames(B4_fluxes)])))
ggplot2::ggplot(plt_df, ggplot2::aes(x=uptake, y=biomass_diff)) +
  ggplot2::geom_point() + ggplot2::geom_smooth(method='lm', se=FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::xlab('Oxygen Uptake (mmol/gDW/h)') +
  ggplot2::ylab('Absolute Biomass Difference (mmol/gDW/h)')


# 3. No Glutamine:
B82_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/medium_changes/B8.2_and_B9_biomass.csv',
                         row.names=1, check.names=FALSE)
B82_fluxes = 2 * ((1 / (1 + exp(-B82_fluxes_raw))) - 0.5)
plot_df = data.frame(values=c(as.numeric(pfba_normalBlood_biomass['MAR13082', colnames(B82_fluxes)]),
                              as.numeric(B82_fluxes['MAR13082', colnames(B82_fluxes)])),
                     cell_types=rep(metadata[colnames(B82_fluxes), 'cell_type'], 2),
                     state=rep(metadata[colnames(B82_fluxes), 'state'], 2),
                     medium=factor(c(rep('Normal', dim(B82_fluxes)[2]),
                                     rep('No Glutamine', dim(B82_fluxes)[2])),
                                   levels=c('Normal', 'No Glutamine')),
                     paired=c(1:174, 1:174))

ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey',
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')

# 3.1. DNA production
plot_df = data.frame(values=c(as.numeric(pfba_normalBlood_biomass['MAR07160', colnames(B82_fluxes)]),
                              as.numeric(B82_fluxes['MAR07160', colnames(B82_fluxes)])),
                     cell_types=rep(metadata[colnames(B82_fluxes), 'cell_type'], 2),
                     state=rep(metadata[colnames(B82_fluxes), 'state'], 2),
                     medium=factor(c(rep('Normal', dim(B82_fluxes)[2]),
                                     rep('No Glutamine', dim(B82_fluxes)[2])),
                                   levels=c('Normal', 'No Glutamine')),
                     paired=c(1:174, 1:174))
ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey',
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')

# 3.2. RNA production
plot_df = data.frame(values=c(as.numeric(colSums(pfba_normalBlood_biomass[c('MAR07161', 'MAR07162'),
                                                            colnames(B82_fluxes)])),
                              as.numeric(colSums(B82_fluxes[c('MAR07161', 'MAR07162'),
                                                            colnames(B82_fluxes)]))),
                     cell_types=rep(metadata[colnames(B82_fluxes), 'cell_type'], 2),
                     state=rep(metadata[colnames(B82_fluxes), 'state'], 2),
                     medium=factor(c(rep('Normal', dim(B82_fluxes)[2]),
                                     rep('No Glutamine', dim(B82_fluxes)[2])),
                                   levels=c('Normal', 'No Glutamine')),
                     paired=c(1:174, 1:174))
ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey',
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')


# 4. No Nucleotides:
# 4.1. Biomass
B17_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/medium_changes/B17_biomass.csv',
                          row.names=1, check.names=FALSE)
B17_fluxes = 2 * ((1 / (1 + exp(-B17_fluxes_raw))) - 0.5)
plot_df = data.frame(values=c(as.numeric(pfba_normalBlood_biomass['MAR13082', colnames(B17_fluxes)]),
                              as.numeric(B17_fluxes['MAR13082', colnames(B17_fluxes)])),
                     cell_types=rep(metadata[colnames(B17_fluxes),'cell_type'], 2),
                     state=rep(metadata[colnames(B17_fluxes),'state'], 2),
                     medium=factor(c(rep('Normal', dim(B17_fluxes)[2]),
                                     rep('No Nucleotides', dim(B17_fluxes)[2])),
                                   levels=c('Normal', 'No Nucleotides')),
                     paired=c(1:174, 1:174))

ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey',
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')

# 4.2. Production of nucleotides
plot_df = data.frame(values=c(as.numeric(pfba_normalBlood_biomass['MAR07160', colnames(B17_fluxes)]),
                              as.numeric(B17_fluxes['MAR07160', colnames(B17_fluxes)])),
                     cell_types=rep(metadata[colnames(B17_fluxes), 'cell_type'], 2),
                     state=rep(metadata[colnames(B17_fluxes), 'state'], 2),
                     medium=factor(c(rep('Normal', dim(B17_fluxes)[2]),
                                     rep('No Nucleotides', dim(B17_fluxes)[2])),
                                   levels=c('Normal', 'No Nucleotides')),
                     paired=c(1:174, 1:174))

ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey',
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')


# 5. No Glucose:
B85_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/medium_changes/B8.5_biomass.csv',
                          row.names=1, check.names=FALSE)
B85_fluxes = 2 * ((1 / (1 + exp(-B85_fluxes_raw))) - 0.5)
plot_df = data.frame(values=c(as.numeric(pfba_normalBlood_biomass['MAR13082', colnames(B85_fluxes)]),
                              as.numeric(B85_fluxes['MAR13082', colnames(B85_fluxes)])),
                     cell_types=rep(metadata[colnames(B85_fluxes), 'cell_type'], 2),
                     state=rep(metadata[colnames(B85_fluxes), 'state'], 2),
                     medium=factor(c(rep('Normal', dim(B85_fluxes)[2]),
                                     rep('No Glucose', dim(B85_fluxes)[2])),
                                   levels=c('Normal', 'No Glucose')),
                     paired=c(1:174, 1:174))

ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey',
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')


# 6. No chloride:
chloride_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/medium_changes/chloride_biomass.csv',
                           row.names=1, check.names=FALSE)
chloride_fluxes = 2 * ((1 / (1 + exp(-chloride_fluxes_raw))) - 0.5)
plot_df = data.frame(values=c(as.numeric(pfba_normalBlood_biomass['MAR13082', colnames(chloride_fluxes)]),
                              as.numeric(chloride_fluxes['MAR13082', colnames(chloride_fluxes)])),
                     cell_types=rep(metadata[colnames(chloride_fluxes), 'cell_type'], 2),
                     state=rep(metadata[colnames(chloride_fluxes), 'state'], 2),
                     medium=factor(c(rep('Normal', dim(chloride_fluxes)[2]),
                                     rep('No Glucose', dim(chloride_fluxes)[2])),
                                   levels=c('Normal', 'No Glucose')),
                     paired=c(1:174, 1:174))

ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey',
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')


# 6. No Serine + Glycine:
B11_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/medium_changes/B11.csv',
                          row.names=1, check.names=FALSE)
B11_fluxes = 2 * ((1 / (1 + exp(-B11_fluxes_raw))) - 0.5)
plot_df = data.frame(values=c(as.numeric(fba_normalBlood['MAR13082', rownames(metadata)]),
                              as.numeric(B11_fluxes['MAR13082', rownames(metadata)])),
                     cell_types=rep(metadata[,'cell_type'], 2),
                     state=rep(metadata[,'state'], 2),
                     medium=factor(c(rep('Normal', dim(metadata)[1]),
                                     rep('No Serine and Glycine', dim(metadata)[1])),
                                   levels=c('Normal', 'No Serine and Glycine')),
                     paired=c(1:196, 1:196))

ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey',
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')


# 7. No arginine:
B121_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/medium_changes/B12.1.csv',
                          row.names=1, check.names=FALSE)
B121_fluxes = 2 * ((1 / (1 + exp(-B121_fluxes_raw))) - 0.5)
plot_df = data.frame(values=c(as.numeric(fba_normalBlood['MAR13082', rownames(metadata)]),
                              as.numeric(B121_fluxes['MAR13082', rownames(metadata)])),
                     cell_types=rep(metadata[,'cell_type'], 2),
                     state=rep(metadata[,'state'], 2),
                     medium=factor(c(rep('Normal', dim(metadata)[1]),
                                     rep('No Arginine', dim(metadata)[1])),
                                   levels=c('Normal', 'No Arginine')),
                     paired=c(1:196, 1:196))

ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey',
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')


# 8. No cysteine:
B122_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/medium_changes/B12.2.csv',
                           row.names=1, check.names=FALSE)
B122_fluxes = 2 * ((1 / (1 + exp(-B122_fluxes_raw))) - 0.5)
plot_df = data.frame(values=c(as.numeric(fba_normalBlood['MAR13082', rownames(metadata)]),
                              as.numeric(B122_fluxes['MAR13082', rownames(metadata)])),
                     cell_types=rep(metadata[,'cell_type'], 2),
                     state=rep(metadata[,'state'], 2),
                     medium=factor(c(rep('Normal', dim(metadata)[1]),
                                     rep('No Cysteine', dim(metadata)[1])),
                                   levels=c('Normal', 'No Cysteine')),
                     paired=c(1:196, 1:196))

ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey',
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')


# 9. No cysteine and arginine:
B123_fluxes_raw = read.csv('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/medium_changes/B12.3.csv',
                           row.names=1, check.names=FALSE)
B123_fluxes = 2 * ((1 / (1 + exp(-B123_fluxes_raw))) - 0.5)
plot_df = data.frame(values=c(as.numeric(fba_normalBlood['MAR13082', rownames(metadata)]),
                              as.numeric(B123_fluxes['MAR13082', rownames(metadata)])),
                     cell_types=rep(metadata[,'cell_type'], 2),
                     state=rep(metadata[,'state'], 2),
                     medium=factor(c(rep('Normal', dim(metadata)[1]),
                                     rep('No Arginine and Cysteine', dim(metadata)[1])),
                                   levels=c('Normal', 'No Arginine and Cysteine')),
                     paired=c(1:196, 1:196))

ggplot2::ggplot(plot_df, ggplot2::aes(x=medium, y=values, colour=cell_types)) +
  ggplot2::geom_boxplot(outlier.shape=NA) +
  ggplot2::geom_line(ggplot2::aes(group=paired), linetype=2, colour='grey',
                     position = ggplot2::position_dodge(0.15)) +
  ggplot2::geom_point(ggplot2::aes(colour=cell_types, group=paired),
                      position = ggplot2::position_dodge(0.2)) +
  ggplot2::facet_wrap(ggplot2::vars(cell_types), scales='free') +
  ggplot2::theme_light() + 
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1, vjust=1, size=8), legend.position='none') +
  ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::xlab('') +
  ggplot2::ylab('Flux (mmol/gDW/h)')




# ----------------------
# --- PATHWAY SCORES ---
# ----------------------

# 0. Load annotations (only considering pathways with more than 5 reactions):
HumanGEM_subsystems_rxns = jsonlite::read_json('./GENERAL/utility_data/subsystems_reactions_mapping.json',
                                               simplifyVector = TRUE)
for(path in names(HumanGEM_subsystems_rxns)){
  rxns = HumanGEM_subsystems_rxns[[path]][HumanGEM_subsystems_rxns[[path]]%in%rownames(fba_normalBlood)]
  if(length(rxns)>5) HumanGEM_subsystems_rxns[[path]] = rxns
  else HumanGEM_subsystems_rxns[[path]] = NULL
}

# 1. Calculate the pathway activity score:
n_models = dim(metadata)[1]
n_pathways = length(HumanGEM_subsystems_rxns)

# 1.1. Create matrix where results are stored:
pathway_scores = matrix(rep(0, n_models*n_pathways), nrow=n_pathways)
dimnames(pathway_scores) = list(names(HumanGEM_subsystems_rxns), rownames(metadata))

# 1.2. Populate matrix:
for(pathway in rownames(pathway_scores)){
  message(pathway)
  for(model in colnames(pathway_scores)){
    path_rxns = HumanGEM_subsystems_rxns[[pathway]]
    path_fluxes = abs(fba_normalBlood[path_rxns, model])
    
    n_rxns = length(path_fluxes)
    ratio_active_rxns = sum(path_fluxes != 0) / n_rxns
    sum_fluxes = sum(path_fluxes)
    path_score = (sum_fluxes / n_rxns) * ratio_active_rxns
    
    pathway_scores[pathway, model] = path_score
  }
}

# 1.3. Heatmap with all pathways:
ht = ComplexHeatmap::Heatmap(pathway_scores, #col = circlize::colorRamp2(c(0, 50, 100), c('blue', 'white', 'red')),
                             name='Score',
                             column_split = metadata[colnames(pathway_scores),]$cell_type, column_title = NULL,
                             row_names_max_width = ComplexHeatmap::max_text_width(names(paths_use), gp = grid::gpar(fontsize = 12)),
                             show_column_names = FALSE, row_names_gp = grid::gpar(fontsize=6),
                             top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata[colnames(pathway_scores),]$cell_type,
                                                                                state=metadata[colnames(pathway_scores),]$state,
                                                                                cms=metadata[colnames(pathway_scores),]$CMS,
                                                                                col=list(cell_type=ct_colors, state=state_colors,
                                                                                         cms=cms_colors),
                                                                                annotation_legend_param = list(cell_type=list(title='Cell Type'),
                                                                                                               state=list(title='State'),
                                                                                                               cms=list(title='CMS'))),
                             heatmap_legend_param = list(direction = "horizontal"))
ComplexHeatmap::draw(ht, merge_legend = TRUE, heatmap_legend_side = 'right', show_heatmap_legend=TRUE)

# 1.3. Heatmap with all pathways:
paths_use = sort(rowMeans(pathway_scores), decreasing = T)[1:20]
ht = ComplexHeatmap::Heatmap(pathway_scores[names(paths_use),], #col = circlize::colorRamp2(c(0, 50, 100), c('blue', 'white', 'red')),
                             name='Score',
                             column_split = metadata[colnames(pathway_scores),]$cell_type, column_title = NULL,
                             row_names_max_width = ComplexHeatmap::max_text_width(names(paths_use), gp = grid::gpar(fontsize = 12)),
                             show_column_names = FALSE, row_names_gp = grid::gpar(fontsize=8),
                             top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadata[colnames(pathway_scores),]$cell_type,
                                                                                state=metadata[colnames(pathway_scores),]$state,
                                                                                cms=metadata[colnames(pathway_scores),]$CMS,
                                                                                col=list(cell_type=ct_colors, state=state_colors,
                                                                                         cms=cms_colors),
                                                                                annotation_legend_param = list(cell_type=list(title='Cell Type'),
                                                                                                               state=list(title='State'),
                                                                                                               cms=list(title='CMS'))),
                             heatmap_legend_param = list(direction = "horizontal"))
ComplexHeatmap::draw(ht, merge_legend = TRUE, heatmap_legend_side = 'bottom', show_heatmap_legend=TRUE)


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
  for(path in rownames(pathway_scores)){
    fc = gtools::foldchange(mean(pathway_scores[path, tum_cts]), mean(pathway_scores[path, norm_cts]))
    if(!is.na(fc) & abs(fc) > 1.5){
      paths_to_test = c(paths_to_test, path)
      fcs_paths = c(fcs_paths, fc)
    }
  }
  if(length(paths_to_test) != 0){
    p_val = sapply(paths_to_test,
                   FUN = function(x){
                     return(wilcox.test(as.numeric(pathway_scores[x, cts_vector]) ~ groups_vector)$p.value)
                   })
    pval_adjust = p.adjust(p_val, method='fdr')
    res_ct = data.frame(pathway=paths_to_test, pval_adjust, pval=p_val, fold_change=fcs_paths, cell_type=ct)
    diff_pathways = rbind(diff_pathways, res_ct)
  }
}
write.csv(diff_pathways, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/FBA/normalBlood_normalVStumour.csv')

# 1.6. Heatmap for each cell-type of the above results:
hetamaps_normVStum = list()
for(ct in unique(diff_pathways$cell_type)){
  ct_diffs = na.omit(diff_pathways$pathway[diff_pathways$pval_adjust < 0.05 & diff_pathways$cell_type==ct])
  ct_models = rownames(metadata)[metadata$cell_type==ct]
  if(length(ct_diffs) == 1){
    mat = t(as.data.frame(pathway_scores[ct_diffs, ct_models]))
    rownames(mat) = ct_diffs
  } 
  else mat = pathway_scores[ct_diffs, ct_models]
  hetamaps_normVStum[[ct]] = ComplexHeatmap::Heatmap(mat,
                                                     #col = circlize::colorRamp2(c(0, 50, 100), c('blue', 'white', 'red')),
                                                     name='Score',
                                                     #column_km = 3,
                                                     row_names_max_width = ComplexHeatmap::max_text_width(rownames(pathway_scores[ct_diffs, ]), 
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


ggplot2::ggplot(data.frame(glycolysis=pathway_scores['Glycolysis / Gluconeogenesis',],
                           oxphos=pathway_scores['Oxidative phosphorylation',],
                           cell_type=metadata[colnames(pathway_scores), 'cell_type'],
                           cms=metadata[colnames(pathway_scores), 'CMS']),
                ggplot2::aes(x=oxphos, y=glycolysis, colour=cell_type)) +
  ggplot2::geom_point(size=3) + ggplot2::scale_colour_manual(values=ct_colors) +
  ggplot2::theme_light()






# ----------------------------
# --- Fluxes similarity ---
# ----------------------------

distance_matrix = as.matrix(dist(t(pfba_normalBlood)))
metadatap = metadata[colnames(pfba_normalBlood),]

ht = ComplexHeatmap::Heatmap(distance_matrix,
                             col = circlize::colorRamp2(c(0, 35, 70), c('blue', 'white', 'red')),
                             name='Euclidean Distance',
                             column_split = metadatap$cell_type, column_title = NULL,
                             row_split = metadatap$cell_type, row_title = NULL,
                             show_column_names = FALSE, show_row_names = FALSE,
                             top_annotation = ComplexHeatmap::HeatmapAnnotation(cell_type=metadatap$cell_type,
                                                                                state=metadatap$state,
                                                                                cms=metadatap$CMS,
                                                                                col=list(cell_type=ct_colors, state=state_colors, cms=cms_colors),
                                                                                annotation_legend_param = list(cell_type=list(title='Cell Type'),
                                                                                                               state=list(title='State'),
                                                                                                               cms=list(title='CMS')),
                                                                                annotation_label = c('Cell Type', 'State', 'CMS')),
                             left_annotation = ComplexHeatmap::rowAnnotation(cell_type=metadatap$cell_type,
                                                                             state=metadatap$state,
                                                                             cms=metadatap$CMS,
                                                                             col=list(cell_type=ct_colors, state=state_colors, cms=cms_colors),
                                                                             annotation_legend_param = list(cell_type=list(title='Cell Type'),
                                                                                                            state=list(title='State'),
                                                                                                            cms=list(title='CMS')),
                                                                             annotation_label = c('Cell Type', 'State', 'CMS')),
                             heatmap_legend_param = list(direction = "horizontal"))
ComplexHeatmap::draw(ht, merge_legend = TRUE, heatmap_legend_side = 'right')#, show_heatmap_legend=FALSE)

