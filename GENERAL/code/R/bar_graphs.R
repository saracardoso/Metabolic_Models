important_rxns = jsonlite::read_json('./GENERAL/utility_data/important_reactions_Tcells.json',
                                     simplifyVector=T)

# -------
# - FAO -
# -------
fadh2 = c('MAR03522', 'MAR03527', 'MAR03531', 'MAR03275', 'MAR03280', 'MAR03284', 'MAR03294', 'MAR03107',
          'MAR03111', 'MAR03115', 'MAR03121', 'MAR03128', 'MAR03135', 'MAR03142', 'MAR03149', 'MAR03156',
          'MAR03202', 'MAR03194', 'MAR03190', 'MAR03186', 'MAR03182', 'MAR03198', 'MAR03178', 'MAR03174',
          'MAR03170', 'MAR03240', 'MAR03244', 'MAR03250', 'MAR03258', 'MAR03218', 'MAR03222', 'MAR03226',
          'MAR03234')

nadhh = c('MAR03524', 'MAR03529', 'MAR03533', 'MAR03278', 'MAR03282', 'MAR03286', 'MAR03292', 'MAR03304',
          'MAR03309', 'MAR03314', 'MAR03319', 'MAR03109', 'MAR03113', 'MAR03117', 'MAR03123', 'MAR03130',
          'MAR03137', 'MAR03144', 'MAR03151', 'MAR03158', 'MAR03204', 'MAR03488', 'MAR03503', 'MAR03509',
          'MAR03514', 'MAR01220', 'MAR01225', 'MAR01178', 'MAR01183', 'MAR03242', 'MAR03246', 'MAR03254',
          'MAR03200', 'MAR03196', 'MAR03192', 'MAR03188', 'MAR03184', 'MAR03180', 'MAR03176', 'MAR03172',
          'MAR03262', 'MAR03358', 'MAR03362', 'MAR03220', 'MAR03224', 'MAR03228', 'MAR03232', 'MAR03236',
          'MAR03328', 'MAR03332', 'MAR03336', 'MAR03340', 'MAR03344', 'MAR03348', 'MAR03352', 'MAR03372',
          'MAR03321')

acetylcoa = c('MAR03525', 'MAR03534', 'MAR03279', 'MAR03283', 'MAR03287', 'MAR03293', 'MAR03305',
              'MAR03310', 'MAR03315', 'MAR03320', 'MAR03160', 'MAR03110', 'MAR03114', 'MAR03118',
              'MAR03125', 'MAR03132', 'MAR03139', 'MAR03146', 'MAR03153', 'MAR03205', 'MAR03201',
              'MAR03197', 'MAR03193', 'MAR03189', 'MAR03185', 'MAR03181', 'MAR03177', 'MAR03173',
              'MAR03510', 'MAR01216', 'MAR01226', 'MAR01174', 'MAR01174', 'MAR01179', 'MAR01184',
              'MAR03243', 'MAR03247', 'MAR03256', 'MAR03264', 'MAR03359', 'MAR03363', 'MAR03221',
              'MAR03225', 'MAR03229', 'MAR03230', 'MAR03233', 'MAR03237', 'MAR03329', 'MAR03333',
              'MAR03337', 'MAR03341', 'MAR03345', 'MAR03349', 'MAR03353', 'MAR03373', 'MAR03368',
              'MAR03367')

all_fadh2 = c('MAR00965', 'MAR00970', 'MAR00970', 'MAR00923')

all_fnah = c('MAR03422', 'MAR03428', 'MAR03430', 'MAR09719', 'MAR00818', 'MAR00922', 'MAR09719',
             'MAR03427', 'MAR03426', 'MAR00789', 'MAR03421', 'MAR00788', 'MAR03409', 'MAR03431',
             'MAR03413', 'MAR03424', 'MAR03423', 'MAR03416', 'MAR03411', 'MAR03432', 'MAR03429',
             'MAR03433', 'MAR03414', 'MAR03406', 'MAR03398', 'MAR03408', 'MAR03425', 'MAR03407',
             'MAR03396', 'MAR03397')

other_fao = c('MAR03523', 'MAR03526', 'MAR03528', 'MAR03530', 'MAR03532', 'MAR03277', 'MAR03281',
              'MAR03285', 'MAR03288', 'MAR03290', 'MAR03296', 'MAR03298', 'MAR03301', 'MAR03302',
              'MAR03306', 'MAR03307', 'MAR03311', 'MAR03312', 'MAR03316', 'MAR03317', 'MAR03321',
              'MAR03322', 'MAR03323', 'MAR03108', 'MAR03112', 'MAR03116', 'MAR03122', 'MAR03129',
              'MAR03136', 'MAR03143', 'MAR03150', 'MAR03157', 'MAR03171', 'MAR03175', 'MAR03179',
              'MAR03183', 'MAR03187', 'MAR03191', 'MAR03195', 'MAR03199', 'MAR03203', 'MAR03480',
              'MAR03482', 'MAR03484', 'MAR03486', 'MAR03489', 'MAR03491', 'MAR03493', 'MAR03498',
              'MAR03501', 'MAR03505', 'MAR03506', 'MAR03508', 'MAR03511', 'MAR03512', 'MAR03513',
              'MAR03515', 'MAR01217', 'MAR01218', 'MAR01219', 'MAR01223', 'MAR01222', 'MAR01224',
              'MAR03241', 'MAR03245', 'MAR03252', 'MAR03260', 'MAR03272', 'MAR03356', 'MAR03357',
              'MAR03360', 'MAR03361', 'MAR03219', 'MAR03223', 'MAR03227', 'MAR03231', 'MAR03235',
              'MAR03239', 'MAR03326', 'MAR03327', 'MAR03330', 'MAR03331', 'MAR03334', 'MAR03335',
              'MAR03338', 'MAR03339', 'MAR03342', 'MAR03343', 'MAR03346', 'MAR03347', 'MAR03350',
              'MAR03351', 'MAR03355', 'MAR03375', 'MAR03370', 'MAR03369', 'MAR03365', 'MAR03364')

fao_bar_graph = function(fba_fluxes, sample){
  total_all_fadh2_flux = sum(fba_fluxes[all_fadh2, sample])
  total_all_fnah_flux = sum(fba_fluxes[all_fnah, sample])
  
  total_fadh2_flux = sum(fba_fluxes[fadh2, sample]) + total_all_fnah_flux
  total_nadhh_flux = sum(fba_fluxes[nadhh, sample]) + total_all_fadh2_flux + total_all_fnah_flux
  xx=fba_fluxes[acetylcoa, sample]
  xx[xx<0] = 0
  total_acetylcoa_flux = sum(xx) + total_all_fadh2_flux + total_all_fnah_flux
  plot_df = data.frame(fluxes=c(total_fadh2_flux, total_nadhh_flux, total_acetylcoa_flux),
                       production=c('FADH2', 'NADH', 'Acetyl-CoA'))
  plot_df$production = factor(plot_df$production, levels=c('FADH2', 'NADH', 'Acetyl-CoA'))
  ggplot2::ggplot(plot_df, ggplot2::aes(production, fluxes)) +
    ggplot2::geom_bar(stat='identity', fill='steelblue', color='black') +
    ggplot2::xlab('') + ggplot2::ylab('Flux (mmol/gDW/h)') +
    ggplot2::theme_minimal()
}



# --------------------------------------------------------------------
# - FADH2 vs NADH between TCA vs glycolysis vs glutaminolysis vs FAO -
# --------------------------------------------------------------------

nadh_fadh_bar_graph = function(fba_fluxes, sample){
  total_all_fadh2_flux = sum(fba_fluxes[all_fadh2, sample])
  total_all_fnah_flux = sum(fba_fluxes[all_fnah, sample])
  fadh2_FAO = sum(fba_fluxes[fadh2, sample]) + total_all_fnah_flux
  nadhh_FAO = sum(fba_fluxes[nadhh, sample]) + total_all_fadh2_flux + total_all_fnah_flux
  
  fadh2_TCA = fba_fluxes['MAR08743', sample]
  if(fadh2_TCA > 0) fadh2_TCA = 0
  else fadh2_TCA = abs(fadh2_TCA)
  xx = fba_fluxes[c('MAR04588', 'MAR04141'), sample]
  xx[xx>0] = 0
  xx = abs(xx)
  nadh_TCA = sum(c(xx, fba_fluxes[c('MAR03957', 'MAR05297'), sample]))
  
  xx = fba_fluxes[c('MAR04373', 'MAR04388'), sample]
  xx[xx>0] = 0
  xx = abs(xx)
  nadh_glycolysis = sum(xx)
  
  nadh_glutaminolysis = fba_fluxes['MAR03802', sample]
  if(nadh_glutaminolysis>0) nadh_glutaminolysis = 0
  else nadh_glutaminolysis = abs(nadh_glutaminolysis)
  
  plot_df = data.frame(fluxes = c(fadh2_FAO, nadhh_FAO, fadh2_TCA, nadh_TCA, nadh_glycolysis,
                                  nadh_glutaminolysis),
                       production = c('FADH2', 'NADH', 'FADH2', rep('NADH', 3)),
                       pathway = c(rep('FAO', 2), rep('TCA', 2), 'Glycolysis', 'Glutaminolysis'))
  plot_df$production = factor(plot_df$production, levels=c('FADH2', 'NADH'))
  plot_df$pathway = factor(plot_df$pathway, levels=c('FAO', 'TCA', 'Glycolysis', 'Glutaminolysis'))
  ggplot2::ggplot(plot_df, ggplot2::aes(x=production, y=fluxes, fill=pathway)) +
    ggplot2::geom_bar(stat='identity', position=ggplot2::position_dodge(), color='black') +
    ggplot2::xlab('') + ggplot2::ylab('Flux (mmol/gDW/h)') +
    ggplot2::theme_minimal()
}



# --------------------------
# - FAO vs FA uptake vs FAS-
# --------------------------

fao_fau_fas_bar_graph = function(fba_fluxes, sample){
  # Total flux trough FAO:
  total_all_fadh2_flux = sum(fba_fluxes[all_fadh2, sample])
  total_all_fnah_flux = sum(fba_fluxes[all_fnah, sample])
  fadh2_flux = sum(fba_fluxes[fadh2, sample])
  nadhh_flux = sum(fba_fluxes[nadhh, sample])
  xx=fba_fluxes[acetylcoa, sample]
  xx[xx<0] = 0
  acetylcoa_flux = sum(xx)
  xx=fba_fluxes[other_fao, sample]
  names(xx) = other_fao
  xx[xx[c('MAR03298', 'MAR03323', 'MAR03480', 'MAR03482', 'MAR03491', 'MAR03355', 'MAR03375')]>0] = 0
  xx[xx[c('MAR03298', 'MAR03323', 'MAR03480', 'MAR03482', 'MAR03491', 'MAR03355', 'MAR03375')]>0] =
    abs(xx[xx[c('MAR03298', 'MAR03323', 'MAR03480', 'MAR03482', 'MAR03491', 'MAR03355', 'MAR03375')]>0])
  xx[xx[c('MAR03489', 'MAR03493', 'MAR01222')]<0] = 0
  other_fao_flux = sum(xx)
  total_fao = total_all_fadh2_flux + total_all_fnah_flux + fadh2_flux + nadhh_flux + acetylcoa_flux + other_fao_flux
  # Total flux of FA uptake:
  total_fau = sum(na.omit(fba_fluxes[important_rxns$uptakes$`fatty acids`, sample]))
  # Total flux through FAS:
  total_fas = sum(na.omit(fba_fluxes[important_rxns$subsystems$`Fatty acid biosynthesis`, sample]))
  # Plot bar graph
  plot_df = data.frame(fluxes=c(total_fao, total_fau, total_fas),
                       path=c('FAO', 'FA Uptake', 'FAS'))
  plot_df$path = factor(plot_df$path, levels=c('FA Uptake', 'FAO', 'FAS'))
  ggplot2::ggplot(plot_df, ggplot2::aes(path, fluxes)) +
    ggplot2::geom_bar(stat='identity', fill='steelblue', color='black') +
    ggplot2::xlab('') + ggplot2::ylab('Flux (mmol/gDW/h)') +
    ggplot2::theme_minimal()
}



