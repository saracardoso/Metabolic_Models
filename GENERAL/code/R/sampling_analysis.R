# Density plots of a reaction flux separated by sample state and cell-type:
density_rxn_StateCT = function(data, metadata, rxn,
                               colour_by=NULL, cts_order=NULL,
                               state_order=NULL){
  # Create data.frame for plotting:
  samplings_names = rownames(data)
  df = cbind(data[samplings_names, rxn], metadata[samplings_names, ])
  colnames(df)[1] = rxn
  if(is.null(cts_order)) df$cell_type = factor(df$cell_type,
                                               levels=unique(df$cell_type))
  else df$cell_type = factor(df$cell_type, levels=cts_order)
  if(is.null(state_order)) df$state = factor(df$state, levels=unique(df$state))
  else df$state = factor(df$state, levels=state_order)
  if(!is.null(colour_by)) df[,colour_by] = factor(df[,colour_by])
  
  # Create plot:
  if(is.null(colour_by)) plt = ggplot2::ggplot(df, ggplot2::aes_string(rxn))
  else plt = ggplot2::ggplot(df, ggplot2::aes_string(rxn, colour=colour_by))
  plt = plt + ggplot2::geom_density() +
    ggplot2::facet_grid(rows=ggplot2::vars(state),
                        cols=ggplot2::vars(cell_type))
  
  return(plt)
}

# Density plots of several reactions whose fluxes are summed, separated by
# sample state and cell-type:
density_sumrxns_StateCT = function(data, metadata, rxns, rxns_rep,
                               colour_by=NULL, cts_order=NULL,
                               state_order=NULL, abs=FALSE){
  # Create data.frame for plotting:
  samplings_names = rownames(data)
  if(abs) aggregated_fluxes = Matrix::rowMeans(abs(data[samplings_names, rxns]))
  else aggregated_fluxes = Matrix::rowMeans(data[samplings_names, rxns])
  df = cbind(aggregated_fluxes, metadata[samplings_names, ])
  colnames(df)[1] = rxns_rep
  if(is.null(cts_order)) df$cell_type = factor(df$cell_type,
                                               levels=unique(df$cell_type))
  else df$cell_type = factor(df$cell_type, levels=cts_order)
  if(is.null(state_order)) df$state = factor(df$state, levels=unique(df$state))
  else df$state = factor(df$state, levels=state_order)
  
  # Create plot:
  if(is.null(colour_by)) plt = ggplot2::ggplot(df, ggplot2::aes_string(rxns_rep))
  else plt = ggplot2::ggplot(df, ggplot2::aes_string(rxns_rep, colour=colour_by))
  plt = plt + ggplot2::geom_density() +
    ggplot2::facet_grid(rows=ggplot2::vars(state),
                        cols=ggplot2::vars(cell_type))
  
  return(plt)
}





















# ---------------------------------------------------------------------------


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
  plt = ggplot2::ggplot(umap_result, ggplot2::aes(UMAP_1, UMAP_2))
  if(!is.null(shape_meta)) 
    plt = plt + ggplot2::geom_point(ggplot2::aes_string(colour=colour_meta, shape=shape_meta))
  else plt = plt + ggplot2::geom_point(ggplot2::aes_string(colour=colour_meta))
  return(plt)
}

# Differential expression analysis




