# Histogram plot for the reaction's fluxes from a sampling result
histogram_sampling = function(sampling_result, reaction){
  plot = ggplot2::ggplot(sampling_result, ggplot2::aes(reaction)) +
    ggplot2::geom_histogram(binwidth=0.01) +
    ggplot2::xlab('Flux (mmol gDW-1 hr-1)') + ggplot2::ylab('Nº Samples')
  return(plot)
}

# Frequency lines plot of multiple reactions' fluxes from a sampling result
frequency_sampling = function(sampling_result, reactions){
  fluxes = c()
  reactions = c()
  for(rxn in reactions){
    fluxes = c(fluxes, sampling_result[,rxn])
    reactions = c(reactions, rep(rxn, dim(sampling_result)[1]))
  }
  df = data.frame(fluxes, reactions)
  plot = ggplot2::ggplot(df, ggplot2::aes(fluxes, colour=reactions)) +
    ggplot2::geom_freqpoly(binwidth=0.01) +
    ggplot2::xlab('Flux (mmol gDW-1 hr-1)') + ggplot2::ylab('Nº Samples')
  return(plot)
}

# Frequency lines plot of multiple reactions' fluxes for each sampling result
frequency_sampling = function(sampling_results, reactions){
  fluxes = c()
  reactions = c()
  sampling_evaluation = c()
  for(samp_eval in names(sampling_results)){
    for(rxn in reactions){
      fluxes = c(fluxes, sampling_results[[samp_eval]][,rxn])
      reactions = c(reactions, rep(rxn, dim(sampling_result[[samp_eval]])[1]))
      sampling_evaluation = c(sampling_evaluation, rep(samp_eval, dim(sampling_result[[samp_eval]])[1]))
    }
  }
  
  df = data.frame(fluxes, reactions, sampling_evaluation)
  plot = ggplot2::ggplot(df, ggplot2::aes(fluxes, colour=reactions)) + ggplot2::facet_wrap(ggplot2::vars(sampling_evaluation)) +
    ggplot2::geom_freqpoly(binwidth=0.01) +
    ggplot2::xlab('Flux (mmol gDW-1 hr-1)') + ggplot2::ylab('Nº Samples')
  return(plot)
}

# UMAP 
umap_sampling = function(sampling_results){
  data_for_umap = c()
  for(samp_eval in names(sampling_results)){
    data_add = sampling_results[[samp_eval]]
    rownames(data_add) = paste(samp_eval, rownames(data_add), sep='_')
    data_for_umap = rbind(data_for_umap, data_add)
  }
  res = umap::umap(data_for_umap)
  return(res)
}

# metadata should be a table with info about cell-type model, state, sample, individual, evaluation sampling.
# umap_result: cbind(umap result, metadata dataframe) -> metadata variables must be factors!
# colour_meta: metadata variable name to colour the samples by
# shape_meta: metadata variable to shape the samples by
plot_umap = function(umap_result, colour_meta, shape_meta=NULL){ 
  # numeric data is: UMAP_1: umap_result$layout[,1]; UMAP_2: umap_result$layout[,2]
  df_umap = data.frame(UMAP_1 = umap_result$layout[,1], UMAP_2 = umap_result$layout[,2])
  df_umap = cbind(df_umap, metadata)
  plt = ggplot2::ggplot(df_umap, ggplot2::aes(UMAP_1, UMAP_2))
  if(!is.null(shape_meta)) 
    plt = plt + ggplot2::geom_point(ggplot2::aes_string(colour=colour_meta, shape=shape_meta))
  else plt = plt + ggplot2::geom_point(ggplot2::aes_string(colour=colour_meta))
  return(plt)
}

# Differential expression analysis




