CRCatlas_file = './0Data/scRNAseq/CRC_atlas/atlas/CRCatlas_PH.h5Seurat'



# ---
# - Read CRC atlas file
# ---

CRCatlas = SeuratDisk::LoadH5Seurat(CRCatlas_file)





# ---
# - Know what samples to use
# ---

# Samples must have a total of at least 1 000 cells.
# For tumours with normal matched mucosa, all samples from an individual must have at least 1 000 cells to consider all samples.
# For healthy donors, samples from each individual will be summed. As such, the sum of all samples should have more than 1 000 cells

annot_4_cells = unique(CRCatlas$Annotation_4[!CRCatlas$Annotation_2%in%c('Normal Epithelial cells', 'IgA+ Plasma cells', 'IgG+ Plasma cells',
                                                                         'Unknown')])
indivs_to_use = as.data.frame(matrix(rep('', length(unique(CRCatlas$Sample)) * 40), ncol=40,
                                     dimnames=list(row=unique(CRCatlas$Sample),
                                                   column=c('Individual', 'Sample.Source', 'Normal Epithelial cells', annot_4_cells,
                                                            'Plasma cells', 'USE'))))
for(samp in unique(CRCatlas$Sample)){
  message(samp)
  indivs_to_use[samp, c('Individual', 'Sample.Source')] = unique(CRCatlas[[]][CRCatlas$Sample==samp, c('Individual', 'Sample.Source')])
  annot_2 = table(CRCatlas[[]][CRCatlas$Sample==samp, 'Annotation_2'])['Normal Epithelial cells']
  annot_2[is.na(annot_2)] = 0
  annot_4 = table(CRCatlas[[]][CRCatlas$Sample==samp, 'Annotation_4'])[annot_4_cells]
  annot_4[is.na(annot_4)] = 0
  
  indivs_to_use[samp, c('Normal Epithelial cells')] = annot_2
  indivs_to_use[samp, annot_4_cells] = annot_4
  
  plasma = table(CRCatlas[[]][CRCatlas$Sample==samp, 'Annotation_2'])[c('IgA+ Plasma cells', 'IgG+ Plasma cells')]
  plasma[is.na(plasma)] = 0
  indivs_to_use[samp, 'Plasma cells'] = sum(plasma)
  
  if(indivs_to_use[samp, 'Sample.Source']!='Healthy Donors'){
    if(sum(as.numeric(indivs_to_use[samp, 3:39])) < 1000) indivs_to_use[samp, 'USE'] = FALSE
  }
}

for(indiv in unique(indivs_to_use$Individual[indivs_to_use$Sample.Source=='Healthy Donors'])){
  samps = rownames(indivs_to_use)[indivs_to_use$Individual==indiv]
  if(sum(as.numeric(unlist(indivs_to_use[samps, 3:39]))) < 1000) indivs_to_use[samps, 'USE'] = FALSE
}

indivs_not_to_use = c('36', '37', '38', 'KUL28', 'KUL30', 'KUL31', 'SMC02', 'SMC03', 'SMC05', 'SMC09', 'N13')





# ---
# - Aggregate gene counts by cell-types and samples, and get metadata
# ---

# A cell-type will only be kept as part of the sample if it is represented by at least 5 cells.
# A gene will only be kept if is present in at least one of the cell-types for that sample.
# IgA+ Plasma cells and IgG+ Plasma cells will be joined into Plasma cells
# Epithelial cells will only be separated into Normal and Tumour cells
# All other cells will have annotation_4

CRCatlas[['models_annotation']] = CRCatlas$Annotation_4
CRCatlas$models_annotation[CRCatlas$Annotation_2=='Normal Epithelial cells'] = 'Normal Epithelial cells'
CRCatlas$models_annotation[CRCatlas$Annotation_4%in%c('IgA+ Plasma cells', 'IgG+ Plasma cells')] = 'Plasma cells'
expression_data_samples = list()
indivs_meta = c()
samps_meta = c()
samp_source_meta = c()
ncells_meta = c()
for(indiv in unique(CRCatlas$Individual)){
  message(indiv)
  if(indiv%in%indivs_not_to_use) next
  indiv_samples = unique(CRCatlas$Sample[CRCatlas$Individual==indiv])
  if(length(indiv_samples)==4){ # Healthy Donors
    indiv_cells = rownames(CRCatlas[[]])[CRCatlas$Individual==indiv]
    aggregated_exp = Seurat::AggregateExpression(subset(CRCatlas, cells=indiv_cells), assays='RNA', slot='counts',
                                                 group.by = 'models_annotation')$RNA
    
    nCells_per_ct = table(CRCatlas[[]][indiv_cells, 'models_annotation'])
    nCells_per_ct = nCells_per_ct[names(nCells_per_ct)!='Unknown']
    final = aggregated_exp[rowSums(aggregated_exp)!=0, names(nCells_per_ct)[nCells_per_ct>5]]
    expression_data_samples[[indiv]] = final
    
    ncells_samp = rep(0, length(unique(CRCatlas$models_annotation)))
    names(ncells_samp) = unique(CRCatlas$models_annotation)
    ncells_samp[names(nCells_per_ct)] = nCells_per_ct
    ncells_meta = rbind(ncells_meta, ncells_samp)
    samp_source_meta = c(samp_source_meta, unique(CRCatlas[[]][CRCatlas$Individual==indiv, 'Sample.Source']))
    indivs_meta = c(indivs_meta, indiv)
    samps_meta = c(samps_meta, indiv)
    
    invisible(gc())
  }
  else{
    expression_data_samples[[indiv]] = list()
    for(samp in indiv_samples){
      message(paste('-', samp))
      samp_cells = rownames(CRCatlas[[]])[CRCatlas$Sample==samp]
      aggregated_exp = Seurat::AggregateExpression(subset(CRCatlas, cells=samp_cells), assays='RNA', slot='counts',
                                                   group.by = 'models_annotation')$RNA
      
      nCells_per_ct = table(CRCatlas[[]][samp_cells, 'models_annotation'])
      nCells_per_ct = nCells_per_ct[names(nCells_per_ct)!='Unknown']
      final = aggregated_exp[rowSums(aggregated_exp)!=0, names(nCells_per_ct)[nCells_per_ct>5]]
      expression_data_samples[[indiv]][[samp]] = final
      
      ncells_samp = rep(0, length(unique(CRCatlas$models_annotation)))
      names(ncells_samp) = unique(CRCatlas$models_annotation)
      ncells_samp[names(nCells_per_ct)] = nCells_per_ct
      ncells_meta = rbind(ncells_meta, ncells_samp)
      samp_source_meta = c(samp_source_meta, unique(CRCatlas[[]][CRCatlas$Sample==samp, 'Sample.Source']))
      indivs_meta = c(indivs_meta, indiv)
      samps_meta = c(samps_meta, samp)
      
      invisible(gc())
    }
  }
}
metadata = data.frame(Individual=indivs_meta, Sample.Source=samp_source_meta)
metadata = cbind(metadata, ncells_meta)
rownames(metadata) = samps_meta





# ---
# - Save data and metadata
# ---

# 1. Save data:

for(indiv in names(expression_data_samples)){
  message(indiv)
  if(unique(metadata[metadata$Individual==indiv,'Sample.Source'])=='Healthy Donors'){
    write.csv(expression_data_samples[[indiv]], paste('./0Data/scRNAseq/CRC_atlas/expression_data/Healthy/', indiv, '.csv', sep=''))
  }
  else{
    if(length(expression_data_samples[[indiv]])>1){ # With normal matched
      dir.create(paste('./0Data/scRNAseq/CRC_atlas/expression_data/CRC/normal_matched', indiv, sep='/'))
      for(samp in names(expression_data_samples[[indiv]])){
        message(samp)
        write.csv(expression_data_samples[[indiv]][[samp]],
                  paste('./0Data/scRNAseq/CRC_atlas/expression_data/CRC/normal_matched/', indiv, '/', samp, '.csv', sep=''))
      }
    }
    else{ # Without normal matched
      write.csv(expression_data_samples[[indiv]][[1]],
                paste('./0Data/scRNAseq/CRC_atlas/expression_data/CRC/no_normal_matched/', names(expression_data_samples[[indiv]]),
                      '.csv', sep=''))
    }
  }
}


# 2. Save metadata:
write.csv(metadata, './0Data/scRNAseq/CRC_atlas/expression_data/metadata.csv')


