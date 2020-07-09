#Directories:
base_dir = '~/Documents/PhD'

Metabolic_Models_Repo = paste(base_dir, 'Metabolic_Models', sep='/')
Metabolic_Models_Repo_general_code = paste(Metabolic_Models_Repo, 'general_code', sep='/')
Metabolic_Models_Repo_general_utils = paste(Metabolic_Models_Repo_general_code, 'utils', sep='/')

#Files:
recon3DModel_genes_file = paste(Metabolic_Models_Repo_general_utils, 'recon3DModel_genes.csv', sep='/')
recon3DModel_GPR_file = paste(Metabolic_Models_Repo_general_utils, 'recon3DModel_GPR.txt', sep='/')





#############
###PA CALLS##
#############

#-- MAIN FUNCTION





###########################

#---READ DATA FOR PA CALLS

#ensembl_genes    -->   vector with the ensembl gene ids to consider. Must correspond to those that are present in the
#                       generic model that will be used for reconstruction
PA_prepare_data = function(ensembl_genes,
                                 transcriptomics_data_files = list(RNAseq=NULL, Microarray=NULL),
                                 transcriptomics_metadata_files = list(RNAseq=NULL, Microarray=NULL),
                                 proteomics_data_files = list(dataset_name=NULL),
                                 proteomics_metadata_files = list(dataset_name=NULL)){
  if(is.null(transcriptomics_data_files$RNAseq)&is.null(transcriptomics_data_files$Microarray))
    stop("Please give a file path to the RNAseq data (transcriptomics_data_files$RNAseq) and/or
         a file path to the Microarray data (transcriptomics_data_files$Microarray)")
  
  expression_data = list()
  expression_data$Transcriptomics = list()
  expression_data$Proteomics = list()
  
  #Get transcriptomics data and maintain only genes in recond3DModel:
  for(transcriptomics_type in names(transcriptomics_data_files)){
    if(!transcriptomics_type%in%c('RNAseq','Microarray')) stop(paste('Invalid transcriptomics type (',
                                                                     transcriptomics_type,
                                                                     '). Valid types: RNAseq and Microarray',
                                                                     sep=''))
    expression_data$Transcriptomics[[transcriptomics_type]] = list()
    data_df = read.csv(transcriptomics_data_files[[transcriptomics_type]], row.names=1)
    if(transcriptomics_type=='RNAseq') expression_data$Transcriptomics[[transcriptomics_type]]$data = data_df[ensembl_genes,]
    else expression_data$Transcriptomics[[transcriptomics_type]]$data = data_df[na.omit(match(ensembl_genes, rownames(data_df))),]
    if(is.null(transcriptomics_metadata_files[[transcriptomics_type]])) warning(paste('Transcriptomics type (',
                                                                                      transcriptomics_type,
                                                                                      ') does not have a corresponding metadata file path under same name',
                                                                                      sep=''))
    expression_data$Transcriptomics[[transcriptomics_type]]$metadata = read.csv(transcriptomics_metadata_files[[transcriptomics_type]],
                                                                row.names = 1)
  }
  
  #Get proteomics data and maintain only genes in recond3DModel:
  for(proteomics_type in names(proteomics_data_files)){
    #if(!proteomics_type%in%c('protein_copy_numbers', 'protein_concentration')) stop(paste('Invalid proteomics type (',
    #                                                                                      proteomics_type,
    #                                                                                      '). Valid types: protein_copy_numbers and protein_concentration',
    #                                                                                      sep=''))
    expression_data$Proteomics[[proteomics_type]] = list()
    data_df = read.csv(proteomics_data_files[[proteomics_type]], row.names=1)
    expression_data$Proteomics[[proteomics_type]]$data = data_df[na.omit(match(ensembl_genes, rownames(data_df))),]
    if(is.null(transcriptomics_metadata_files[[proteomics_type]])) warning(paste('Proteomics data (',
                                                                                 proteomics_type,
                                                                                 ') does not have a corresponding metadata file path under same name',
                                                                                 sep=''))
    else expression_data$Proteomics[[proteomics_type]]$metadata = read.csv(transcriptomics_metadata_files[[proteomics_type]],
                                                                           row.names = 1)
  }
  
  return(expression_data)
}



#---GENE PA CALLS: TRANSCRIPTOMICS

#gene_mapping   -->   data.frame with columns 'ensembl', 'entrez' and 'symbol'. Corresponds to the genes that are present
#                     in the model to be used as generic model in the reconstruction
PA_transcriptomics_calls = function(transcriptomics_data, gene_mapping){
  
  #PA calls will be stored in this dataset:
  if(!'RNAseq'%in%names(transcriptomics_data)) RNAseq_init = rep(NA,length(gene_mapping$entrez))
  else RNAseq_init = rep(0,length(gene_mapping$entrez))
  if(!'Microarray'%in%names(transcriptomics_data)) Microarray_init = rep(NA,length(gene_mapping$entrez))
  else Microarray_init = rep(0,length(gene_mapping$entrez))
  final_gene_calls = data.frame(RNAseq=RNAseq_init, Microarray=Microarray_init,
                                Transcriptomics=rep(-1,length(gene_mapping$entrez)),
                                row.names=gene_mapping$entrez)
  
  #List where PA results will be stored:
  PA_res = list()
  PA_res$calls_per_sample = list()
  
  #PA calls for each data type:
  for(exp_data in names(transcriptomics_data)){
    
    if(!is.null(transcriptomics_data[[exp_data]]$data)){
      n_genes = dim(transcriptomics_data[[exp_data]]$data)[1]
      n_samples = dim(transcriptomics_data[[exp_data]]$data)[2]
      PA_calls = matrix(rep(0,n_genes*n_samples), nrow=n_genes, ncol=n_samples)
      colnames(PA_calls) = colnames(transcriptomics_data[[exp_data]]$data)
      rownames(PA_calls) = rownames(transcriptomics_data[[exp_data]]$data)
      
      dens = summary(as.vector(as.matrix(transcriptomics_data[[exp_data]]$data)))
      overall_Q1 = dens[2]
      overall_Q3 = dens[5]
      
      for(gene in rownames(transcriptomics_data[[exp_data]]$data)){
        gene_mean = mean(as.numeric(transcriptomics_data[[exp_data]]$data[gene,]))
        t = overall_Q3
        if(gene_mean < overall_Q3) t = max(c(overall_Q1, gene_mean))
        PA_calls[gene, transcriptomics_data[[exp_data]]$data[gene,]>=t] = 1
      }
      
      PA_res$calls_per_sample[[exp_data]] = as.data.frame(matrix(rep(0,length(gene_mapping$entrez)*n_samples),
                                  nrow=length(gene_mapping$entrez), ncol=n_samples))
      colnames(PA_res$calls_per_sample[[exp_data]]) = colnames(transcriptomics_data[[exp_data]]$data)
      rownames(PA_res$calls_per_sample[[exp_data]]) = gene_mapping$entrez
      for(samp in colnames(PA_res$calls_per_sample[[exp_data]])){
        pa_cals_intermediate = PA_calls[,samp]
        names(pa_cals_intermediate) = rownames(PA_calls)
        PA_res$calls_per_sample[[exp_data]][,samp] = pa_cals_intermediate[gene_mapping$ensembl]
      } 
      
      PA_calls_average = rowMeans(PA_calls)
      final_gene_calls[,exp_data] = PA_calls_average[gene_mapping$ensembl]
    }
    else final_gene_calls[,exp_data] = rep(NA, length(gene_mapping$ensembl))
    
  }
  
  #Final PA calls:
  for (gene in 1:dim(final_gene_calls)[1]){ #For each gene
    
    if(sum(!is.na(final_gene_calls[gene,c('RNAseq','Microarray')]))==2){ #If there are PA calls from both data types
      average_gene = mean(as.numeric(final_gene_calls[gene,c('RNAseq','Microarray')])) #Average the score
      if(average_gene<=0.5) final_gene_calls[gene, 'Transcriptomics'] = 0 #If average smaller than 0.5, score 0 (not present)
      else if(average_gene>=0.75) final_gene_calls[gene, 'Transcriptomics'] = 2 #If average greater than 0.75,
                                                                          #score 2 (high confidence for being present)
      else final_gene_calls[gene, 'Transcriptomics'] = 1 #If between 0.5 and 0.75, score 1 (medium confidence for being present)
    }
    else if(sum(!is.na(final_gene_calls[gene,c('RNAseq','Microarray')]))==1){#If there are PA calls only from one data type
      gene_value = na.omit(as.numeric(final_gene_calls[gene,c('RNAseq','Microarray')])) #Get the score from that data type
      if(gene_value<=0.5) final_gene_calls[gene, 'Transcriptomics'] = 0 #If score smaller than 0.5, score 0 (not present)
      else if(gene_value>=0.8) final_gene_calls[gene, 'Transcriptomics'] = 2 #This time, score must be greater than 0.8 to be
                                                                       #scored 2 (high confidence for being present)
      else final_gene_calls[gene, 'Transcriptomics'] = 1 #If between 0.5 and 0.8, score 1 (medium confidence for being present)
    }
    #If is not in either --> it will be -1
  }
  PA_res$calls_per_omics = final_gene_calls
  
  return(PA_res)
}



#---GENE PA CALLS: PROTEOMICS

#gene_mapping   -->   data.frame with columns 'ensembl', 'entrez' and 'symbol'. Corresponds to the genes that are present
#                     in the model to be used as generic model in the reconstruction
PA_proteomics_calls = function(proteomics_data, gene_mapping){
  
  #PA calls will be stored in this dataset:
  final_gene_calls = c()
  for(exp_data in names(proteomics_data)){
    if(!is.null(proteomics_data[[exp_data]]$data)) final_gene_calls = cbind(final_gene_calls, rep(NA,length(gene_mapping$entrez)))
  }
  final_gene_calls = cbind(final_gene_calls, rep(NA,length(gene_mapping$entrez)))
  colnames(final_gene_calls) = c(names(proteomics_data), 'Proteomics')
  row.names(final_gene_calls) = gene_mapping$entrez
  
  #List where PA results will be stored:
  PA_res = list()
  PA_res$calls_per_sample = list()
  
  #PA calls for each data type:
  for(exp_data in names(proteomics_data)){
    
    if(!is.null(proteomics_data[[exp_data]]$data)){
      n_genes = dim(proteomics_data[[exp_data]]$data)[1]
      n_samples = dim(proteomics_data[[exp_data]]$data)[2]
      PA_calls = matrix(rep(0,n_genes*n_samples), nrow=n_genes, ncol=n_samples)
      colnames(PA_calls) = colnames(proteomics_data[[exp_data]]$data)
      rownames(PA_calls) = rownames(proteomics_data[[exp_data]]$data)
      
      dens = summary(as.vector(as.matrix(proteomics_data[[exp_data]]$data)))
      overall_Q1 = dens[2]
      overall_Q3 = dens[5]
      
      for(gene in rownames(proteomics_data[[exp_data]]$data)){
        gene_mean = mean(as.numeric(proteomics_data[[exp_data]]$data[gene,]))
        t = overall_Q3
        if(gene_mean < overall_Q3) t = max(c(overall_Q1, gene_mean))
        PA_calls[gene, proteomics_data[[exp_data]]$data[gene,]>=t] = 1
      }
      
      PA_res$calls_per_sample[[exp_data]] = as.data.frame(matrix(rep(0,length(gene_mapping$entrez)*n_samples),
                                                nrow=length(gene_mapping$entrez), ncol=n_samples))
      colnames(PA_res$calls_per_sample[[exp_data]]) = colnames(proteomics_data[[exp_data]]$data)
      rownames(PA_res$calls_per_sample[[exp_data]]) = gene_mapping$entrez
      for(samp in colnames(PA_res$calls_per_sample[[exp_data]])){
        pa_cals_intermediate = PA_calls[,samp]
        names(pa_cals_intermediate) = rownames(PA_calls)
        PA_res$calls_per_sample[[exp_data]][,samp] = pa_cals_intermediate[gene_mapping$ensembl]
      } 
      
      PA_calls_average = rowMeans(PA_calls)
      final_gene_calls[,exp_data] = PA_calls_average[gene_mapping$ensembl]
    }
    else final_gene_calls[,exp_data] = rep(NA, length(gene_mapping$ensembl))
    
  }
  
  #Final PA calls:
  for (gene in 1:dim(final_gene_calls)[1]){ #For each gene
    average_gene = mean(as.numeric(final_gene_calls[gene,dim(final_gene_calls)[2]-1])) #Average the score
    if(!is.na(average_gene)){
      if(average_gene<0.5) final_gene_calls[gene, 'Proteomics'] = 0 #If score smaller than 0.5, score 0 (not present)
      else if(average_gene>=0.8) final_gene_calls[gene, 'Proteomics'] = 2 #This time, score must be greater than 0.8 to be
      #scored 2 (high confidence for being present)
      else final_gene_calls[gene, 'Proteomics'] = 1 #If between 0.5 and 0.8, score 1 (medium confidence for being present)
    }
    
  }
  PA_res$calls_per_omics = final_gene_calls
  
  return(PA_res)
}


#---FROM GENE TO REACTION CALLS

#RULES        -->   OR: MAX or SUM; AND: MIN.
#gene_calls   -->   data.frame of one column. Row names: entrez gene ids. First column: gene scores.
#gpr_rules    -->   data.frame with two columns. First: reaction IDs from the generic model. Second: GPR rules
PA_from_gene_to_reaction_calls = function(gene_calls, gpr_rules_filtered, or = 'MAX'){
  
  #Get GPR rules for each reaction in the recon3D model:
  ##PUT THIS IN SEPARATE FUNCTION::
  #recon3D_gpr_rules = read.table(recon3DModel_GPR_file, sep='\t')
  #recon3D_gpr_rules_filtered = recon3D_gpr_rules[!recon3D_gpr_rules[,2]=='',] #Remove reactions without GPRs
  #recon3D_gpr_rules_filtered[,2] = gsub('[.][1-9]+', '',recon3D_gpr_rules_filtered[,2]) #Remove 'version' .1, .2, ...
  
  #Reaction calls:
  reaction_calls = c()
  for(reaction in 1:dim(gpr_rules_filtered)[1]){
    gpr_str = gpr_rules_filtered[reaction,2]
    gpr_str = gsub('[(]','',gpr_str) #Remove (), they won't be necessary.
    gpr_str = gsub('[)]','',gpr_str)
    splitted_ors = strsplit(gpr_str,' or ')[[1]] #Split into different elements through or
    for(i in 1:length(splitted_ors)){ #Get minimum value in each of these elements
      genePA_calls = gene_calls[strsplit(splitted_ors[i],' and ')[[1]], 1]
      splitted_ors[i] = min(genePA_calls)
    }
    if(length(na.omit(splitted_ors))==0) reaction_calls = c(reaction_calls, 0)
    else{
      if(or=='MAX') reaction_calls = c(reaction_calls, max(as.numeric(na.omit(splitted_ors)))) #Between the resulting scores of the elements, get
                                                                                               #maximum value (as these are connected through or)
      else if(or=='SUM') reaction_calls = c(reaction_calls, sum(as.numeric(na.omit(splitted_ors)))) #Between the resulting scores of the elements, sum
                                                                                                    #the values (as these are connected through or)
      else stop('Invalid OR rule. or argument must be MAX or SUM.')
    }
  }
  names(reaction_calls) = gpr_rules_filtered[,1]
  
  return(reaction_calls)
}


#-- FINAL REACTION CALLS

PA_reaction_calls = function(transcriptomics_reaction_calls, proteomics_reaction_calls=NULL){
  
}








final_gene_PA_calls = function(transcriptomics_calls, proteomics_calls){
  
  #Results will be stored here:
  final_gene_calls = list()
  final_gene_calls$calls_per_sample = list()
  final_gene_calls$calls_per_sample$Transcriptomics = transcriptomics_calls$calls_per_sample
  if(!is.null(proteomics_calls)) final_gene_calls$calls_per_sample$Proteomics = proteomics_calls$calls_per_sample
  final_gene_calls$calls_per_omics = list()
  final_gene_calls$calls_per_omics$Transcriptomics = transcriptomics_calls$calls_per_omics
  if(!is.null(proteomics_calls)) final_gene_calls$calls_per_omics$Proteomics = proteomics_calls$calls_per_omics
  
  #Get mapping between entrez, ensembl and gene ids:
  recon3DModel_genes = read.csv(recon3DModel_genes_file, header=T)
  
  #Get final gene PA calls:
  final_calls = c()
  
  if(!is.null(proteomics_calls)){
  for(gene in recon3DModel_genes$entrez){
    transcriptomics_score = final_gene_calls$calls_per_omics$Transcriptomics[as.character(gene),'Transcriptomics']
    proteomics_score = final_gene_calls$calls_per_omics$Proteomics[as.character(gene),'Proteomics']
    if(is.na(transcriptomics_score)&is.na(proteomics_score)) gene_score = 0 # If no data on gene
    else if(sum(is.na(c(transcriptomics_score, proteomics_score)))==1){ # If one of the omics does not have data for the gene
      if(na.omit(c(transcriptomics_score, proteomics_score))==2) gene_score = 1
      else gene_score = 0
    }
    else{ # When both omics have data for the gene
      if(transcriptomics_score+proteomics_score <= 1) gene_score = 0
      else if(abs(transcriptomics_score-proteomics_score) == 2){
        if(transcriptomics_score==2 & sum(na.omit(as.numeric(final_gene_calls$calls_per_omics$Proteomics[as.character(gene),])))==0) gene_score = 0
        else if(proteomics_score==2 & sum(na.omit(as.numeric(final_gene_calls$calls_per_omics$Transcriptomics[as.character(gene),])))==0) gene_score = 0
        else gene_score = 1
      }
      else if(transcriptomics_score+proteomics_score>=3) gene_score = 2
      else gene_score = 1
    }
    final_calls = c(final_calls, gene_score)
  }
  final_gene_calls$final_gene_calls = data.frame(Transcriptomics=final_gene_calls$calls_per_omics$Transcriptomics[,'Transcriptomics'],
                                                 Proteomics=final_gene_calls$calls_per_omics$Proteomics[,'Proteomics'],
                                                 Final=final_calls,
                                                 row.names=recon3DModel_genes$entrez)
  }
  else{
    for(gene in recon3DModel_genes$entrez){
      rnaseq_score = final_gene_calls$calls_per_omics$Transcriptomics[as.character(gene),'RNAseq']
      microarray_score = final_gene_calls$calls_per_omics$Transcriptomics[as.character(gene),'Microarray']
      
      if(is.na(rnaseq_score)&is.na(microarray_score)) gene_score = 0 # If no data on gene
      else if(sum(is.na(c(rnaseq_score, microarray_score)))==1){ # If one of the transcriptomics does not have data for the gene
        if(na.omit(c(rnaseq_score, microarray_score))>=0.8) gene_score = 2
        if(na.omit(c(rnaseq_score, microarray_score))<=0.5) gene_score = 0
        else gene_score = 1
      }
      else{
        if(rnaseq_score>=0.75 & microarray_score>=0.75) gene_score = 2
        else if (mean(c(rnaseq_score, microarray_score))<=0.5) gene_score = 0
        else gene_score = 1
      }
      final_calls = c(final_calls, gene_score)
    }
    final_gene_calls$final_gene_calls = data.frame(Transcriptomics=final_gene_calls$calls_per_omics$Transcriptomics[,'Transcriptomics'],
                                                   Final=final_calls,
                                                   row.names=recon3DModel_genes$entrez)
  }
  
  return(final_gene_calls)
}









##################
####GO ANALYSIS###
##################

GO_analysis = function(interest_genes, universe_genes, org_db='org.Hs.eg.db', keytype='ENSEMBL',
                       pvalue_cutoff = 0.05, p_adjust_method = 'fdr', qvalue_cutoff = 0.10){
  
  final_result = list()
  for(ontology in c('BP','CC','MF')){
    final_result[[ontology]] = clusterProfiler::enrichGO(gene=interest_genes, universe=universe_genes, OrgDb=org_db,
                                           keyType=keytype, ont=ontology, readable=T,
                                           pvalueCutoff = pvalue_cutoff, pAdjustMethod = p_adjust_method,
                                           qvalueCutoff = qvalue_cutoff)
  }
  return(final_result)
}

GO_simplify_result = function(GO_result){
  simple_GO_result = list()
  for(ontology in names(GO_result)){
    ont_res = GO_result[[ontology]]
    simple_GO_result[[ontology]] = clusterProfiler::simplify(ont_res, by='p.adjust')
  }
  return(simple_GO_result)
}

#GO_result    -->   A result from running GO_analysis or simplify_GO_result
significant_GOs = function(GO_result, p.adjust_cutoff = 0.05){
  significant_GOs_tables = list()
  for(ontology in names(GO_result)){
    ont_res = GO_result[[ontology]]
    significant_GOs_tables[[ontology]] = ont_res@result[ont_res@result$p.adjust<p.adjust_cutoff,
                                                        c('Description','GeneRatio','BgRatio','p.adjust','qvalue','Count')]
  }
  return(significant_GOs_tables)
}





########################
####CLUSTER ANALYSIS####
########################

cluster_analysis = function(dataset, scale_data=F, color_leafs = 'GSE', label_leafs = 'data_type', plot_title='', leg_pos=c(70,25)){
  if(scale_data) dataset_dist = dist(scale(t(dataset$data)))
  else dataset_dist = dist(t(dataset$data))
  dataset_hclust = hclust(dataset_dist)
  pallete = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
              "#B15928", "#808080", "#000000")
  dendogram_plot_col_modified(dataset, dataset_hclust, color_leafs, label_leafs, colors = pallete,
                              leg.pos=leg_pos, lab.cex=0.9, title = plot_title)
}

color_leaf_modified = function (n, dataset, classes, colors, labels, lab.cex = 1){
  if (is.leaf(n)) {
    a <- attributes(n)
    i <- match(a$label, specmine::get_sample_names(dataset))
    attr(n, "nodePar") <- c(a$nodePar, list(lab.col = colors[classes[i]], 
                                            lab.cex = lab.cex, pch = NA))
    attr(n, "label") <- labels[i]
  }
  n
}

dendogram_plot_col_modified = function (dataset, hc.result, classes.col, colors = NULL, title = "", 
                                        lab.cex = 1, leg.pos = c(1,1), label_samples = NULL, ...){
  classes = dataset$metadata[, classes.col]
  cluster = as.dendrogram(hc.result)
  if (is.null(label_samples)) 
    label_samples = colnames(dataset$data)
  else label_samples = dataset$metadata[, label_samples]
  if (is.null(colors)) 
    colors = 1:length(levels(classes))
  cluster = dendrapply(cluster, color_leaf_modified, dataset, classes, colors,
                       label_samples, lab.cex)
  plot(cluster, main = title, horiz = FALSE, bty='L', ...)
  leg.txt = levels(classes)
  leg.txt = c("Legend:", leg.txt)
  leg.col = c("black", colors)
  if (!is.null(leg.pos))
    par(xpd=TRUE)
  legend(leg.pos[1], leg.pos[2], leg.txt, text.col = leg.col, bty = "n", cex=lab.cex, y.intersp=.5, x.intersp=.2)
}





########################
####HEATMAP ANALYSIS####
########################

heatmap_gene_calls = function(cell_type_PA_genes, plot_tilte='',
                              sample_colors_legend_location=c(1,1), cell_colours_legend_location=c(1,2)){#, preclusterrows=F){
  pallete = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
              "#B15928", "#000000")
  
  #Organize data for heatmap and set sample colors:
  exp_data_present = c(names(cell_type_PA_genes$calls_per_sample$Transcriptomics),
                       names(cell_type_PA_genes$calls_per_sample$Proteomics))
  sample_colors = c()
  i=1
  cell_type_heatmap_data = c()
  for(omics_type in names(cell_type_PA_genes$calls_per_sample)){
    for(type_data in names(cell_type_PA_genes$calls_per_sample[[omics_type]])){
      if(length(cell_type_heatmap_data)==0) cell_type_heatmap_data = cell_type_PA_genes$calls_per_sample[[omics_type]][[type_data]]
      else cell_type_heatmap_data = cbind(cell_type_heatmap_data, cell_type_PA_genes$calls_per_sample[[omics_type]][[type_data]])
      sample_colors = c(sample_colors, rep(pallete[i], ncol(cell_type_PA_genes$calls_per_sample[[omics_type]][[type_data]])))
      i=i+1
    }
  }
  cell_type_heatmap_data = cbind(cell_type_heatmap_data, cell_type_PA_genes$final_gene_calls$Final)
  colnames(cell_type_heatmap_data)[dim(cell_type_heatmap_data)[2]] = 'Final Call'
  rownames(cell_type_heatmap_data) = rownames(cell_type_PA_genes$final_gene_calls)
  cell_type_heatmap_data = as.matrix(cell_type_heatmap_data)
  sample_colors = c(sample_colors, pallete[i])
  
  #Set cells colors:
  col_cells <- c('#ffffff', '#ffc100', '#ff0000', '#808080')
  
  #Order genes according to decreasing presence throughout samples and then by final gene calls:
  #if(preclusterrows){
  #  Rowv <- as.dendrogram(hclust(dist(na.omit(cell_type_heatmap_data))))
  #  rowInd <- order.dendrogram(Rowv)
  #  cell_type_heatmap_data = cell_type_heatmap_data[rowInd,]
  #}
  #cell_type_heatmap_data = cell_type_heatmap_data[order(cell_type_heatmap_data[,'Final Call'], decreasing=T),]
  na_row_sums = rowSums(is.na(cell_type_heatmap_data[,1:dim(cell_type_heatmap_data)[2]-1]))
  cell_type_heatmap_data = cell_type_heatmap_data[order(na_row_sums, decreasing=T),]
  row_sums = rowSums(cell_type_heatmap_data[,1:dim(cell_type_heatmap_data)[2]-1])
  cell_type_heatmap_data = cell_type_heatmap_data[order(row_sums, decreasing=T),]
  cell_type_heatmap_data = cell_type_heatmap_data[order(cell_type_heatmap_data[,'Final Call'], decreasing=T),]
  cell_type_heatmap_data[is.na(cell_type_heatmap_data)] = 3
  
  #Plot heatmap:
  heatmap(cell_type_heatmap_data, Colv=NA, Rowv=NA, scale='none',
          ColSideColors=sample_colors, labRow='', labCol='',
          col=col_cells, cexRow=.2, main=plot_tilte, na.rm=F)
  legend(x=sample_colors_legend_location[1], y=sample_colors_legend_location[2],
         legend=c(exp_data_present, 'Final Call'), fill=pallete[1:i],
         title='Sample colours', cex=.8, bty='n', y.intersp=0.7, x.intersp=0.1, title.adj=0)
  legend(x=cell_colours_legend_location[1], y=cell_colours_legend_location[2],
         legend=c("0         ", "1", "2", "NA"), fill=col_cells,
         title='Cell colours', cex=.8, bty='n', y.intersp=0.7, x.intersp=0.1, title.adj=0)
}





#############
####OTHERS###
#############

reactions_no_GPR = function(){
  #Get GPR rules for each reaction in the recon3D model:
  recon3D_gpr_rules = read.table(recon3DModel_GPR_file, sep='\t')
  reactions_without_GPR = recon3D_gpr_rules[recon3D_gpr_rules[,2]=='',1] #Get reactions without GPRs
  return(reactions_without_GPR)
}






#summary(na.omit(abs(final_gene_calls$RNAseq-final_gene_calls$Microarray)))
#final_gene_calls[na.omit(abs(final_gene_calls$RNAseq-final_gene_calls$Microarray))>.5,]
