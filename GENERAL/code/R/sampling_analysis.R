# Density plots of a reaction flux separated by sample state and cell-type:
boxplot_rxn_StateCT = function(data, metadata, rxn, x,
                               colour_by=NULL, x_order=NULL,
                               colour_order=NULL){
  # Create data.frame for plotting:
  samplings_names = rownames(data)
  df = cbind(data[samplings_names, rxn], metadata[samplings_names, ])
  colnames(df)[1] = rxn
  if(is.null(x_order)) df$cell_type = factor(df$cell_type,
                                               levels=unique(df$cell_type))
  else df$cell_type = factor(df$cell_type, levels=x_order)
  if(is.null(colour_order)) df$state = factor(df$state, levels=unique(df$state))
  else df$state = factor(df$state, levels=colour_order)
  if(!is.null(colour_by)) df[,colour_by] = factor(df[,colour_by])
  
  # Create plot:
  if(is.null(colour_by)) plt = ggplot2::ggplot(df, ggplot2::aes_string(y=rxn, x=x))
  else plt = ggplot2::ggplot(df, ggplot2::aes_string(y=rxn, x=x, colour=colour_by))
  plt = plt + ggplot2::geom_boxplot() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, vjust=.5))#, hjust=1))
  
  return(plt)
}

# Density plots of several reactions whose fluxes are summed, separated by
# sample state and cell-type:
boxplot_sumrxns_StateCT = function(data, metadata, rxns, rxns_rep, x,
                                   colour_by=NULL, x_order=NULL,
                                   colour_order=NULL, abs=FALSE){
  # Create data.frame for plotting:
  samplings_names = rownames(data)
  if(abs) aggregated_fluxes = Matrix::rowMeans(abs(data[samplings_names, rxns]))
  else aggregated_fluxes = Matrix::rowMeans(data[samplings_names, rxns])
  df = cbind(aggregated_fluxes, metadata[samplings_names, ])
  colnames(df)[1] = rxns_rep
  if(is.null(x_order)) df$cell_type = factor(df$cell_type,
                                             levels=unique(df$cell_type))
  else df$cell_type = factor(df$cell_type, levels=x_order)
  if(is.null(colour_order)) df$state = factor(df$state, levels=unique(df$state))
  else df$state = factor(df$state, levels=colour_order)
  if(!is.null(colour_by)) df[,colour_by] = factor(df[,colour_by])
  
  # Create plot:
  if(is.null(colour_by)) plt = ggplot2::ggplot(df, ggplot2::aes_string(y=rxns_rep, x=x))
  else plt = ggplot2::ggplot(df, ggplot2::aes_string(y=rxns_rep, x=x, colour=colour_by))
  plt = plt + ggplot2::geom_boxplot() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, vjust=.5))#, hjust=1))
  
  return(plt)
}

# Density plots of pathway scores separated by sample state
# and cell-type:
density_pathscores_StateCT = function(path_scores_vec, metadata, path_name,
                                      colour_by=NULL, cts_order=NULL,
                                      state_order=NULL){
  # Create data.frame for plotting:
  samplings_names = names(path_scores_vec)
  df = cbind(path_scores_vec[samplings_names], metadata[samplings_names, ])
  colnames(df)[1] = path_name
  if(is.null(cts_order)) df$cell_type = factor(df$cell_type,
                                               levels=unique(df$cell_type))
  else df$cell_type = factor(df$cell_type, levels=cts_order)
  if(is.null(state_order)) df$state = factor(df$state, levels=unique(df$state))
  else df$state = factor(df$state, levels=state_order)
  if(!is.null(colour_by)) df[,colour_by] = factor(df[,colour_by])
  
  # Create plot:
  if(is.null(colour_by)) plt = ggplot2::ggplot(df, ggplot2::aes_string(path_name))
  else plt = ggplot2::ggplot(df, ggplot2::aes_string(path_name, colour=colour_by))
  plt = plt + ggplot2::geom_density() +
    ggplot2::facet_grid(rows=ggplot2::vars(state),
                        cols=ggplot2::vars(cell_type))
  
  return(plt)
}



# Pathway scores:

glycolysis_score = function(data){
  
  glycolysis_rxns = c('MAR04394', 'MAR07747', 'MAR04381', 'MAR04301', 'MAR04379',
                      'MAR04375', 'MAR04391', 'MAR04373', 'MAR04371', 'MAR04372',
                      'MAR04368', 'MAR04370', 'MAR04365', 'MAR04363', 'MAR04378')
  subdata = as.matrix(data[,glycolysis_rxns])
  # glucose -> G6P
  glucose_G6P =rowSums(subdata[,c('MAR04394', 'MAR07747')])
  # G6P -> F6P
  which_pos = subdata[,'MAR04381'] > 0
  subdata[which_pos,'MAR04381'] = 0
  G6P_F6P = abs(subdata[,'MAR04381'])
  # F6P -> F1,6BP
  F6P_F16BP = rowSums(subdata[,c('MAR04301', 'MAR04379')])
  # F1,6BP -> GAP + DHAP
  which_pos = subdata[,'MAR04375'] > 0
  subdata[which_pos,'MAR04375'] = 0
  F16BP_GAPDHAP = abs(subdata[,'MAR04375'])
  # DHAP -> GAP
  which_neg = subdata[,'MAR04391'] < 0
  subdata[which_neg,'MAR04391'] = 0
  DHAP_GAP = subdata[,'MAR04391']
  # GAP -> 1,3BPG
  which_pos = subdata[,'MAR04373'] > 0
  subdata[which_pos,'MAR04373'] = 0
  GAP_13BPG = abs(subdata[,'MAR04373'])
  # 1,3BPG -> 2,3BG
  BPG_23BG = subdata[,'MAR04371']
  # 2,3BG -> 3PG
  BG_3PG = subdata[,'MAR04372']
  # 1,3BPG -> 3PG
  which_neg = subdata[,'MAR04368'] < 0
  subdata[which_neg,'MAR04368'] = 0
  BPG_3PG = rowSums(subdata[,c('MAR04368', 'MAR04370')])
  # 3PG -> 2PD
  which_pos = subdata[,'MAR04365'] > 0
  subdata[which_pos,'MAR04365'] = 0
  PG_2PD = abs(subdata[,'MAR04365'])
  # 2PD -> PEP
  which_neg = subdata[,'MAR04363'] < 0
  subdata[which_neg,'MAR04363'] = 0
  PD_PEP = subdata[,'MAR04363']
  # PEP -> pyruvate
  PEP_pyruvate = subdata[,'MAR04378']
  # Final Score:
  final_scores = rowMeans(cbind(glucose_G6P, G6P_F6P, F6P_F16BP, F16BP_GAPDHAP,
                                DHAP_GAP, GAP_13BPG, BPG_23BG, BG_3PG, BPG_3PG,
                                PG_2PD, PD_PEP, PEP_pyruvate))
  names(final_scores) = rownames(subdata)
  
  return(final_scores)
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




