
###################
##READ OMICS DATA##
###################

#--
#-read_omics_data
#--
#-> Arguments:
# entrez_genes_universe           -->   vector with the entrez gene ids to consider. Must correspond to those that are present
#                                       in the generic model that will be used for reconstruction
# gene_mapping_file               -->   string with the path to the file that contains the differents IDs for at least the genes
#                                       in the generic model. File with 3 columns: entrez, ensembl and symbol
# map_file_sep                    -->   string with character that separates the different cells in the file gene_mapping_file
# map_file_header                 -->   boolean value indicating whether the gene_mapping_file contains a row as header or not
# transcriptomics_data_files      -->   list with two items: RNAseq (string with the path to the RNAseq data file) and
#                                       Microarray (string with the path to the Microarray data file)
# transcriptomics_metadata_files  -->   (optional) list with two items: RNAseq (string with the path to the RNAseq metadata
#                                       file) and Microarray (string with the path to the Microarray metadata file)
# proteomics_data_files           -->   (optional) list with one or more items: each item (named after the dataset's name)
#                                       contains a string with the path to dataset's data file
# proteomics_metadata_files       -->   (optional) list with two items: each item (named after the dataset's name) contains
#                                       a string with the path to dataset's metadata file
#-> Value:
#List with the following structure:
#   -dataset$
#        -Transcriptomics$
#                 -RNAseq$ [if available]
#                     -data
#                     -metadata [if available]
#                 -Microarray$ [if available]
#                     -data
#                     -metadata [if available]
#        -Proteomics$
#                 -<dataset_name>$ [if available]
#                     -data
#                     -metadata [if available]
read_omics_data = function(entrez_genes_universe, gene_mapping_file, map_file_sep=',', map_file_header=TRUE,
                           transcriptomics_data_files = list(RNAseq=NULL, Microarray=NULL),
                           transcriptomics_metadata_files = list(RNAseq=NULL, Microarray=NULL),
                           proteomics_data_files = list(dataset_name=NULL),
                           proteomics_metadata_files = list(dataset_name=NULL)){
  if(is.null(transcriptomics_data_files$RNAseq)&is.null(transcriptomics_data_files$Microarray))
    stop("Please give a file path to the RNAseq data (transcriptomics_data_files$RNAseq) and/or
         a file path to the Microarray data (transcriptomics_data_files$Microarray)")
  
  # Get gene mapping: entrez (mandatory) - ensembl (mandatory) - symbol (optional)
  gene_mapping = read.csv(gene_mapping_file, sep=map_file_sep, header=map_file_header)
  if(length(entrez_genes_universe) > dim(gene_mapping)[1])
    stop('Too many genes given.')
  if(sum(is.na(match(entrez_genes_universe, gene_mapping$entrez)))>0)
    stop('One or more of the genes given are not present in the map_file given.')
  ensembl_genes = gene_mapping$ensembl[match(entrez_genes_universe, gene_mapping$entrez)]
  
  # Data will be stored here:
  expression_data = list()
  expression_data$Transcriptomics = list()
  if(!is.null(proteomics_data_files) | !is.null(proteomics_data_files$dataset_name)) expression_data$Proteomics = list()
  
  #Get transcriptomics data and maintain only genes in recond3DModel:
  for(transcriptomics_type in names(transcriptomics_data_files)){
    if(!transcriptomics_type%in%c('RNAseq','Microarray')) stop(paste('Invalid transcriptomics type (',
                                                                     transcriptomics_type,
                                                                     '). Valid types: RNAseq and Microarray',
                                                                     sep=''))
    expression_data$Transcriptomics[[transcriptomics_type]] = list()
    data_df = read.csv(transcriptomics_data_files[[transcriptomics_type]], row.names=1)
    expression_data$Transcriptomics[[transcriptomics_type]]$data = data_df[na.omit(match(ensembl_genes, rownames(data_df))),]
    if(is.null(transcriptomics_metadata_files[[transcriptomics_type]])) warning(paste('Transcriptomics type (',
                                                                                      transcriptomics_type,
                                                                                      ') does not have a corresponding metadata file path under same name',
                                                                                      sep=''))
    expression_data$Transcriptomics[[transcriptomics_type]]$metadata = read.csv(transcriptomics_metadata_files[[transcriptomics_type]],
                                                                                row.names = 1)
  }
  
  #Get proteomics data and maintain only genes in recond3DModel:
  if(!is.null(proteomics_data_files) | !is.null(proteomics_data_files$dataset_name)){
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
  }
  
  return(expression_data)
}





#############
###PA CALLS##
#############

#-- MAIN FUNCTION

#--
#-PA_reaction_calls
#--
#-> Arguments:
# omics_data                      -->   dataset retrieved from function read_omics_data
# gene_mapping_file               -->   string with the path to the file that contains the differents IDs for at least the genes
#                                       in the generic model. File with 3 columns: entrez, ensembl and symbol
# gpr_file                        -->   string with the path to the file that contains the GPR rules to the reactions
#                                       from the generic model. File with two columns: reaction IDs and GPR rules in the
#                                       form of: ( geneA and geneB ) or ( geneC and geneD ), for example.
# entrez_genes_universe           -->   vector with the entrez gene ids to consider. Must correspond to those that are present
#                                       in the generic model that will be used for reconstruction
# map_file_sep, gpr_file_sep      -->   string with character that separates the different cells in the files
#                                       gene_mapping_file and gpr_file, respectively
# map_file_header,gpr_file_header -->   boolean value indicating whether the gene_mapping_file and gpr_file, respectively,
#                                       contain a row as header or not
# or                              -->   function to use regarding the OR rule. Either 'MAX' or 'SUM'. Defaults to 'MAX'.
#-> Value:
#List with the following structure:
#   -dataset$
#        -calls_per_sample$
#                 -Transcriptomics$
#                     -RNAseq [if available]
#                     -Microarray [if available]
#                 -Proteomics$ [if available]
#                     -<dataset_name>
#        -gene_scores
#        -reaction_calls
PA_reaction_calls = function(omics_data, gene_mapping_file, gpr_file, entrez_genes_universe,
                             or = 'MAX',
                             map_file_sep=',', map_file_header=TRUE,
                             gpr_file_sep=' ', gpr_file_header=FALSE){
  
  # Get GPR rules for reactions in the generic model:
  gpr_rules = read.table(gpr_file, sep=gpr_file_sep, header=gpr_file_header)
  #gpr_rules[,2] = gsub('[.][1-9]+', '',gpr_rules[,2]) #Remove 'version' .1, .2, ...
  
  # Get gene mapping: entrez (mandatory) - ensembl (mandatory) - symbol (optional)
  gene_mapping = read.csv(gene_mapping_file, sep=map_file_sep, header=map_file_header)
  if(length(entrez_genes_universe) > dim(gene_mapping)[1])
    stop('Too many genes given.')
  if(sum(is.na(match(entrez_genes_universe, gene_mapping$entrez)))>0)
    stop('One or more of the genes given are not present in the map_file given.')
  gene_mapping = gene_mapping[match(entrez_genes_universe, gene_mapping$entrez),]
  
  # Gene calls:
  if(is.null(omics_data$Proteomics)){
    #Transcriptomics calls:
    transcriptomics_calls = PA_transcriptomics_calls(omics_data$Transcriptomics, gene_mapping, only_trans=T)
    #Get named vector:
    transcriptomics_gene_calls = transcriptomics_calls$calls_per_omics[,'Transcriptomics']
    names(transcriptomics_gene_calls) = rownames(transcriptomics_calls$calls_per_omics)
    #Get final gene calls:
    gene_calls = PA_gene_calls(transcriptomics_gene_calls)
  }
  else{
    #Transcriptomics calls:
    transcriptomics_calls = PA_transcriptomics_calls(omics_data$Transcriptomics, gene_mapping, only_trans=F)
    #Get named vector:
    transcriptomics_gene_calls = transcriptomics_calls$calls_per_omics[,'Transcriptomics']
    names(transcriptomics_gene_calls) = rownames(transcriptomics_calls$calls_per_omics)
    #Get transcriptomics presence:
    transcriptomics_presence = c()
    for(gene in rownames(transcriptomics_calls$calls_per_omics))
      transcriptomics_presence = c(transcriptomics_presence,
                                   mean(na.omit(as.numeric(transcriptomics_calls$calls_per_omics[gene,1:(dim(transcriptomics_calls$calls_per_omics)[2]-1)]))))
    names(transcriptomics_presence) = rownames(transcriptomics_calls$calls_per_omics)
    
    #Proteomics calls:
    proteomics_calls = PA_proteomics_calls(omics_data$Proteomics, gene_mapping)
    #Get named vector:
    proteomics_gene_calls = proteomics_calls$calls_per_omics[,'Proteomics']
    names(proteomics_gene_calls) = rownames(proteomics_calls$calls_per_omics)
    #Get proteomics presence:
    if(dim(proteomics_calls$calls_per_omics)[2]==2) # Only one proteomics dataset
      proteomics_presence = proteomics_calls$calls_per_omics[,1]
    else{ # More than one proteomics dataset
      proteomics_presence = c()
      for(gene in rownames(proteomics_calls$calls_per_omics))
        proteomics_presence = c(proteomics_presence,
                                mean(na.omit(as.numeric(proteomics_calls$calls_per_omics[gene,1:(dim(proteomics_calls$calls_per_omics)[2]-1)]))))
    }
    
    #Get final gene calls:
    gene_calls = PA_gene_calls(transcriptomics_gene_calls, transcriptomics_presence=transcriptomics_presence,
                               proteomics_gene_calls=proteomics_gene_calls, proteomics_presence=proteomics_presence)
  }
  
  # Reaction calls:
  reaction_calls = PA_from_gene_to_reaction_calls(gene_calls, gpr_rules, or = or)
  
  # Reaction calls for each algorithm:
  algorithm_calls = PA_algorithms_reaction_calls(reaction_calls)
  
  # Store Final Results:
  results = list()
  results$calls_per_sample = list()
  results$calls_per_sample$Transcriptomics = transcriptomics_calls$calls_per_sample
  if(!is.null(omics_data$Proteomics))
    results$calls_per_sample$Proteomics = proteomics_calls$calls_per_sample
  results$gene_scores = gene_calls
  results$reaction_calls = algorithm_calls
  
  return(results)
}




###########################


#---GENE PA CALLS: TRANSCRIPTOMICS

PA_transcriptomics_calls = function(transcriptomics_data, gene_mapping, only_trans=T){
  
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
  
  # Final PA calls:
  
  if(!only_trans){ #If there will be proteomics scores too
               #Here, scores will be -1, 0, 1, 2
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
  }
  else{ #If there will not be proteomics scores to use, transcriptomics scores will be more tight
        #Here, scores will be -1, 0, 1, 2, 3
    for (gene in 1:dim(final_gene_calls)[1]){ #For each gene
      if(sum(!is.na(final_gene_calls[gene,c('RNAseq','Microarray')]))==2){ #If there are PA calls from both data types
        sum_gene = sum(as.numeric(final_gene_calls[gene,c('RNAseq','Microarray')])) #Summ the scores of both types
        average_gene = mean(as.numeric(final_gene_calls[gene,c('RNAseq','Microarray')]))  #Average the score
        if(average_gene<=0.5) final_gene_calls[gene, 'Transcriptomics'] = 0 #If average smaller than 0.5, score 0
        else if(sum_gene>=1.5) final_gene_calls[gene, 'Transcriptomics'] = 3 #If both types greater than 0.75,
                                                                            #score 3
        else final_gene_calls[gene, 'Transcriptomics'] = 1 #Other cases, score 1
      }
      else if(sum(!is.na(final_gene_calls[gene,c('RNAseq','Microarray')]))==1){#If there are PA calls only from one data type
        gene_value = na.omit(as.numeric(final_gene_calls[gene,c('RNAseq','Microarray')])) #Get the score from that data type
        if(gene_value<=0.5) final_gene_calls[gene, 'Transcriptomics'] = 0 #If score smaller than 0.5, score 0
        else if(gene_value>=0.8) final_gene_calls[gene, 'Transcriptomics'] = 2 #This time, score must be greater than 0.8 to be
                                                                              #scored 2
        else final_gene_calls[gene, 'Transcriptomics'] = 1 #If between 0.5 and 0.8, score 1
      }
      #If is not in either --> it will be -1
    }
  }
  
  PA_res$calls_per_omics = final_gene_calls
  
  return(PA_res)
}



#---GENE PA CALLS: PROTEOMICS

PA_proteomics_calls = function(proteomics_data, gene_mapping){
  
  #PA calls will be stored in this dataset:
  final_gene_calls = c()
  for(exp_data in names(proteomics_data)){
    if(!is.null(proteomics_data[[exp_data]]$data)) final_gene_calls = cbind(final_gene_calls, rep(NA,length(gene_mapping$entrez)))
  }
  final_gene_calls = cbind(final_gene_calls, rep(-1,length(gene_mapping$entrez)))
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
  if(length(names(proteomics_data))==1){ # If only one dataset available
    for (gene in 1:dim(final_gene_calls)[1]){ #For each gene
      score_gene = as.numeric(final_gene_calls[gene,dim(final_gene_calls)[2]-1])
      if(!is.na(score_gene)){
        if(score_gene<0.5) final_gene_calls[gene, 'Proteomics'] = 0 #If score smaller than 0.5, score 0 (not present)
        else if(score_gene>=0.8) final_gene_calls[gene, 'Proteomics'] = 2 #This time, score must be greater than 0.8 to be
        #scored 2 (high confidence for being present)
        else final_gene_calls[gene, 'Proteomics'] = 1 #If between 0.5 and 0.8, score 1 (medium confidence for being present)
      }
    #If no data on gene: score is -1
    }
  }
  else{ # If more than one dataset available
    for (gene in 1:dim(final_gene_calls)[1]){ #For each gene
      if(sum(!is.na(final_gene_calls[gene,1:dim(final_gene_calls)[2]-1]))>=2){ #If there are PA calls at least two data types
        average_gene = mean(as.numeric(na.omit(final_gene_calls[gene,1:dim(final_gene_calls)[2]-1]))) #Average the score
        if(average_gene<=0.5) final_gene_calls[gene, 'Transcriptomics'] = 0 #If average smaller than 0.5, score 0 (not present)
        else if(average_gene>=0.75) final_gene_calls[gene, 'Transcriptomics'] = 2 #If average greater than 0.75,
                                                                                  #score 2 (high confidence for being present)
        else final_gene_calls[gene, 'Transcriptomics'] = 1 #If between 0.5 and 0.75, score 1 (medium confidence for being present)
      }
      else if(sum(!is.na(final_gene_calls[gene,1:dim(final_gene_calls)[2]-1]))==1){#If there are PA calls only from one data type
        gene_value = na.omit(as.numeric(final_gene_calls[gene,c('RNAseq','Microarray')])) #Get the score from that data type
        if(gene_value<=0.5) final_gene_calls[gene, 'Transcriptomics'] = 0 #If score smaller than 0.5, score 0 (not present)
        else if(gene_value>=0.8) final_gene_calls[gene, 'Transcriptomics'] = 2 #This time, score must be greater than 0.8 to be
        #scored 2 (high confidence for being present)
        else final_gene_calls[gene, 'Transcriptomics'] = 1 #If between 0.5 and 0.8, score 1 (medium confidence for being present)
      }
      #If is not in either --> it will be -1
    }
  }
  PA_res$calls_per_omics = final_gene_calls
  
  return(PA_res)
}


#---FINAL GENE CALLS

PA_gene_calls = function(transcriptomics_gene_calls, transcriptomics_presence=NULL, 
                         proteomics_gene_calls=NULL, proteomics_presence = NULL){
  
  # IDs between transcriptomics and proteomics (if available) must match:
  if(!is.null(proteomics_gene_calls)){
    if(sum(is.na(match(names(proteomics_gene_calls), names(transcriptomics_gene_calls))))>0)
      stop('Gene IDs do not match entirely between the teo omics.')
    if(is.null(proteomics_presence) | is.null(transcriptomics_presence))
      stop('transcriptomics_presence and proteomics_presence must be given when both proteomics and transcriptomics
           score are given.')
  }
  
  # Results will be stored here:
  final_gene_scores = c()
  
  # In case both omics are available:
  if(!is.null(proteomics_gene_calls)){
    for(gene in names(transcriptomics_gene_calls)){
      both_scores = c(transcriptomics_gene_calls[gene], proteomics_gene_calls[gene])
      scores_sum = sum(both_scores)
      if(transcriptomics_gene_calls[gene]!=-1 & proteomics_gene_calls[gene]!=-1){
        if(scores_sum == 4) final_gene_scores = c(final_gene_scores, 3)
        else if(scores_sum == 3) final_gene_scores = c(final_gene_scores, 2)
        else if(scores_sum <= 1) final_gene_scores = c(final_gene_scores, 0)
        else if(sum(both_scores == c(1,1))==2) final_gene_scores = c(final_gene_scores, 1)
        else{ #When one of the omics is 2 and the other 0
          samps_presence = c(transcriptomics_presence[gene], proteomics_presence[gene])[both_scores==0]
          if(samps_presence==0) final_gene_scores = c(final_gene_scores, 0)
          else final_gene_scores = c(final_gene_scores, 1)
        }
      }
      else if(transcriptomics_gene_calls[gene]!=-1 & proteomics_gene_calls[gene]!=-1) final_gene_scores = c(final_gene_scores, -1)
      else{
        if(scores_sum == 1) final_gene_scores = c(final_gene_scores, 1)
        else final_gene_scores = c(final_gene_scores, 0)
      }
    }
  }
  else final_gene_scores = transcriptomics_gene_calls
  
  names(final_gene_scores) = names(transcriptomics_gene_calls)
  return(final_gene_scores)
}


#---FROM GENE TO REACTION CALLS

#RULES        -->   OR: MAX or SUM; AND: MIN.
#gene_calls   -->   named vector. Names: entrez gene ids. Values: gene scores.
#gpr_rules    -->   data.frame with two columns. First: reaction IDs from the generic model. Second: GPR rules
PA_from_gene_to_reaction_calls = function(gene_calls, gpr_rules, or = 'MAX'){
  
  #Reaction calls:
  reaction_calls = c()
  for(reaction in 1:dim(gpr_rules)[1]){
    gpr_str = gpr_rules[reaction,2]
    gpr_str = gsub('[(] ','',gpr_str) #Remove (), they won't be necessary.
    gpr_str = gsub(' [)]','',gpr_str)
    splitted_ors = strsplit(gpr_str,' or ')[[1]] #Split into different elements through or
    if(length(splitted_ors)==0) reaction_calls = c(reaction_calls, NA) # Reaction with no gene associated
    else{
      for(i in 1:length(splitted_ors)){ #Get minimum value in each of these elements
        genePA_calls = gene_calls[strsplit(splitted_ors[i],' and ')[[1]]]
        splitted_ors[i] = min(genePA_calls)
      }
      if(or=='MAX')
        reaction_calls = c(reaction_calls, max(as.numeric(splitted_ors))) #Between the resulting scores of the elements, get
                                                                          #maximum value (as these are connected through or)
      else if(or=='SUM'){
        splitted_ors = splitted_ors[splitted_ors!=-1] #Take out genes with no experimental data to make the sum
        if(length(splitted_ors)==0) #But if there are no genes with experimental data, give score -1
          reaction_calls = c(reaction_calls -1)
        else
          reaction_calls = c(reaction_calls, sum(as.numeric(splitted_ors))) #Between the resulting scores of the elements, sum
                                                                            #the values (as these are connected through or)
      }
      else stop('Invalid OR rule. or argument must be MAX or SUM.')
    }
    
  }
  names(reaction_calls) = gpr_rules[,1]
  
  return(reaction_calls)
}


#-- FINAL REACTION CALLS

PA_algorithms_reaction_calls = function(reaction_calls){
  
  # Get the ids for the reactions.
  reaction_ids = names(reaction_calls)
  
  # Scores will be calculated for the following algorithms:
  algos = c('FastCore','GIMME','IMAT','CORDA','tINIT')
  
  #Initiate data.frame where results will be stored:
  algo_reaction_calls = data.frame(original=reaction_calls,
                              FastCore=rep('-', length(reaction_ids)),
                              GIMME=rep(-1, length(reaction_ids)),
                              IMAT=rep(-1, length(reaction_ids)),
                              CORDA=rep('-', length(reaction_ids)),
                              tINIT=rep(-2, length(reaction_ids)),
                              row.names=reaction_ids)
  
  # Calculate Reaction Scores for each algorithm:
  for(reaction in reaction_ids){
    if(!is.na(reaction_calls[reaction])){
      if(reaction_calls[reaction] == 3)
        algo_reaction_calls[reaction,c('FastCore','GIMME','IMAT','CORDA','tINIT')] = c('Present', 2, 2, 'High', 20)
      else if(reaction_calls[reaction] == 2)
        algo_reaction_calls[reaction,c('FastCore','GIMME','IMAT','CORDA','tINIT')] = c('Present', 2, 2, 'High', 15)
      else if(reaction_calls[reaction] == 1)
        algo_reaction_calls[reaction,c('FastCore','GIMME','IMAT','CORDA','tINIT')] = c('Present', 1, 1, 'Medium', 10)
      else if(reaction_calls[reaction] == 0)
        algo_reaction_calls[reaction,c('FastCore','GIMME','IMAT','CORDA','tINIT')] = c('-', 0, 0, 'Negative', -8)
    }
  }
  
  return(algo_reaction_calls)
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
