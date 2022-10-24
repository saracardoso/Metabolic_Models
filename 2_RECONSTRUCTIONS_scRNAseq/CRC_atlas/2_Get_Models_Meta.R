
indiv_meta = read.csv('./0Data/scRNAseq/CRC_atlas/expression_data/metadata.csv', row.names=1)

indivs = c()
samples = c()
cell_types = c()
states = c()
n_rxnss = c()
models_names = c()
for(indiv in c('31', '32', '33', '35', 'KUL01', 'KUL19', 'KUL21', 'SMC01', 'SMC04', 'SMC06', 'SMC07',
               'SMC08', 'SMC10')){
  dir_indiv = paste('./2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched', indiv, sep='/')
  n_rxns_files = grep('3_biomassAfter', list.files(dir_indiv, recursive=FALSE), value=TRUE)
  gapfill_file = '2_nReactionsGapFill.csv'
  
  samps = gsub('.csv', '', gsub('3_biomassAfter_', '', n_rxns_files))
  
  n_reactions = read.csv(paste(dir_indiv, '/', '2_nReactionsGapFill.csv', sep=''), row.names=1, check.names=FALSE)
  for(samp in samps){
    biomassAfter = read.csv(paste(dir_indiv, '/', '3_biomassAfter_', samp, '.csv', sep=''), row.names=1)
    cts = rownames(biomassAfter)[biomassAfter$Blood_SMDB!=0]
    n_rxns = n_reactions[cts, samp]
    st = indiv_meta[samp, 'Sample.Source']
    if(st!='Normal Matched') st = 'Tumour'
    
    indivs = c(indivs, rep(indiv, length(cts)))
    samples = c(samples, rep(samp, length(cts)))
    cell_types = c(cell_types, cts)
    states = c(states, rep(st, length(cts)))
    n_rxnss = c(n_rxnss, n_rxns)
    models_names = c(models_names, paste(indiv, samp, cts, sep='_'))
  }
}
final_meta = data.frame(individual=indivs, sample=samples, cell_type=cell_types, state=states,
                        n_reactions=n_rxnss, row.names = models_names)
final_meta = final_meta[final_meta$cell_type!='IL22+ CD4 Tcells',]
# Get CMS types:
tumour_meta = read.csv('./0Data/scRNAseq/CRC_atlas/atlas/Epithelial_tumour.csv', row.names=1)[, c('Sample', 'Annotation_2')]
tumour_meta = unique(tumour_meta)
final_meta$CMS = 'Normal Matched'
for(samp_idx in 1:dim(tumour_meta)[1])
  final_meta[final_meta$sample==tumour_meta$Sample[samp_idx], 'CMS'] = tumour_meta$Annotation_2[samp_idx]
write.csv(final_meta, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/NormalMatched/metadata.csv')




