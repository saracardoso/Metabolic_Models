# --------------
# - Glycolysis -
# --------------
glycolysis_rxns = data.frame(id=c('MAR04394', 'MAR07747', 'MAR04381', 'MAR04301', 'MAR04379', 'MAR04375', 'MAR04391', 'MAR04373', 'MAR04371', 'MAR04372', 'MAR04368', 'MAR04370', 'MAR04365', 'MAR04363', 'MAR04358', 'MAR04926', 'MAR04137', 'MAR04388'),
                              label=c('MAR04394', 'MAR07747', 'MAR04381', 'MAR04301', 'MAR04379', 'MAR04375', 'MAR04391', 'MAR04373', 'MAR04371', 'MAR04372', 'MAR04368', 'MAR04370', 'MAR04365', 'MAR04363', 'MAR04358', 'MAR04926', 'MAR04137', 'MAR04388'),
                              shape='box', shadow='FALSE', color.background=c(rep('lightgrey', 15), 'thistle', rep('lightgrey', 2)),
                              color.border='black', font.color='black', font.size=20)
glycolysis_metabs = data.frame(id=c('MAM01965c', 'MAM01968c', 'MAM01845c', 'MAM01841c', 'MAM01690c', 'MAM01939c', 'MAM00247c', 'MAM00569c', 'MAM00913c', 'MAM00674c', 'MAM02696c', 'MAM02819c', 'MAM02819m', 'MAM01261m', 'MAM02403c'),
                               label=c('glucose', 'glucose-6-phosphate', 'fructose-6-phosphate', 'fructose-1,6-bisphosphate', 'DHAP', 'GAP', '1,3-bisphospho-D-glycerate', '2,3-bisphospho-D-glycerate', '3-phospho-D-glycerate', '2-phospho-D-glycerate', 'PEP', 'pyruvate', 'pyruvate', 'acetyl-CoA', 'lactate'),
                               shape='ellipse', shadow='TRUE', color.background=c(rep('lightblue', 12), rep('gold', 2), 'lightblue'),
                               color.border='black', font.color='black', font.size=20)
glycolysis_mrEdges = data.frame(from=c('MAM01965c', 'MAM01965c', 'MAM01968c', 'MAM01845c', 'MAM01845c', 'MAM01841c', 'MAM01690c', 'MAM01939c', 'MAM00247c', 'MAM00247c', 'MAM00247c', 'MAM00569c', 'MAM00913c', 'MAM00674c', 'MAM02696c', 'MAM02819c', 'MAM02819m', 'MAM02819c'),
                                to=c('MAR04394', 'MAR07747', 'MAR04381', 'MAR04301', 'MAR04379', 'MAR04375', 'MAR04391', 'MAR04373', 'MAR04371', 'MAR04368', 'MAR04370', 'MAR04372', 'MAR04365', 'MAR04363', 'MAR04358', 'MAR04926', 'MAR04137', 'MAR04388'),
                                arrows='to', dashes=FALSE, color='green', width=1)
glycolysis_rmEdges = data.frame(from=c('MAR04394', 'MAR07747', 'MAR04381', 'MAR04301', 'MAR04379', 'MAR04375', 'MAR04375', 'MAR04391', 'MAR04373', 'MAR04371', 'MAR04368', 'MAR04370', 'MAR04372', 'MAR04365', 'MAR04363', 'MAR04358', 'MAR04926', 'MAR04137', 'MAR04388'),
                                to=c('MAM01968c', 'MAM01968c', 'MAM01845c', 'MAM01841c', 'MAM01841c', 'MAM01690c', 'MAM01939c', 'MAM01939c', 'MAM00247c', 'MAM00569c', 'MAM00913c', 'MAM00913c', 'MAM00913c', 'MAM00674c', 'MAM02696c', 'MAM02819c', 'MAM02819m', 'MAM01261m', 'MAM02403c'),
                                arrows='to', dashes=FALSE, color='green', width=1)
glycolysis_edges = rbind(glycolysis_mrEdges, glycolysis_rmEdges)

glycolysis_fluxes_graph = function(rxn_fluxes, main=''){
  for(rxn in glycolysis_rxns$id){
    glycolysis_edges$width[glycolysis_edges$from==rxn | glycolysis_edges$to==rxn] = abs(rxn_fluxes[rxn])*5 + 1
    if(rxn_fluxes[rxn]==0){
      glycolysis_edges$color[glycolysis_edges$from==rxn | glycolysis_edges$to==rxn] = 'grey'
    }
    else if(rxn%in%c('MAR04381', 'MAR04375', 'MAR04373', 'MAR04365') & rxn_fluxes[rxn]>0){
      glycolysis_edges$color[glycolysis_edges$from==rxn | glycolysis_edges$to==rxn] = 'red'
    }
    else if(rxn%in%c('MAR04368', 'MAR04388') & rxn_fluxes[rxn]<0) glycolysis_edges$color[glycolysis_edges$from==rxn | glycolysis_edges$to==rxn] = 'red'
}
visNetwork::visNetwork(rbind(glycolysis_metabs, glycolysis_rxns),
                       glycolysis_edges, main=main)
}



# -------
# - TCA -
# -------
tca_rxns = data.frame(id=c('MAR04145', 'MAR04458', 'MAR04589', 'MAR04588', 'MAR03957', 'MAR03958', 'MAR04112', 'MAR05297', 'MAR03787', 'MAR04147', 'MAR04152', 'MAR04652', 'MAR08743', 'MAR04410', 'MAR04141', 'MAR04209', 'MAR06411', 'MAR06413', 'MAR06414'),
                      label=c('MAR04145', 'MAR04458', 'MAR04589', 'MAR04588', 'MAR03957', 'MAR03958', 'MAR04142', 'MAR05297', 'MAR03787', 'MAR04147', 'MAR04152', 'MAR04652', 'MAR08743', 'MAR04410', 'MAR04141', 'MAR04209', 'MAR06411', 'MAR06413', 'MAR06414'),
                      shape='box', shadow='FALSE', color.background='lightgrey',
                      color.border='black', font.color='black', font.size=20)
tca_metabs = data.frame(id=c('MAM01597m', 'MAM01587m', 'MAM02633m', 'MAM01261m', 'MAM01580m', 'MAM02183m', 'MAM02662m', 'MAM01306m', 'MAM02944m', 'MAM02943m', 'MAM01253m', 'MAM01255m', 'MAM01862m', 'MAM02439m', 'MAM02984m', 'MAM00765m', 'MAM02393m', 'MAM02934m', 'MAM01701m'),
                        label=c('CoA', 'citrate', 'OAA', 'acetyl-CoA', 'cis-aconitate', 'isocitrate', 'oxalosuccinate', 'AKG', 'succinyl-CoA', 'succinate', 'acetoacetate', 'acetoacetyl-CoA', 'fumarate', 'malate', 'thiamin-PP', '3-carboxy-1-hydroxypropyl-ThPP', 'lipoamide', 'S-succinyldihydrolipoamide', 'dihydrolipoamide'),
                        shape='ellipse', shadow='TRUE', color.background='gold',
                        color.border='black', font.color='black', font.size=20)
tca_mrEdges = data.frame(from=c('MAM01261m', 'MAM02633m', 'MAM01587m', 'MAM01580m', 'MAM02183m', 'MAM02662m', 'MAM02183m', 'MAM02183m', 'MAM01306m', 'MAM01597m', 'MAM02944m', 'MAM02944m', 'MAM01253m', 'MAM02944m', 'MAM02943m', 'MAM02943m', 'MAM01862m', 'MAM02439m', 'MAM01306m', 'MAM02984m', 'MAM00765m', 'MAM02393m', 'MAM01306m', 'MAM02393m', 'MAM01597m', 'MAM02934m'),
                         to=c('MAR04145', 'MAR04145', 'MAR04458', 'MAR04589', 'MAR04588', 'MAR04112', 'MAR03957', 'MAR03958', 'MAR05297', 'MAR05297', 'MAR04152', 'MAR04147', 'MAR03787', 'MAR03787', 'MAR08743', 'MAR04652', 'MAR04410', 'MAR04141', 'MAR04209', 'MAR04209', 'MAR06413', 'MAR06413', 'MAR06411', 'MAR06411', 'MAR06414', 'MAR06414'),
                         arrows='to', dashes=FALSE, color='green', width=1)
tca_rmEdges = data.frame(from=c('MAR04145', 'MAR04145', 'MAR04458', 'MAR04589', 'MAR04588', 'MAR04112', 'MAR03957', 'MAR03958', 'MAR05297', 'MAR04152', 'MAR04152', 'MAR04147', 'MAR04147', 'MAR03787', 'MAR03787', 'MAR08743', 'MAR04652', 'MAR04410', 'MAR04141', 'MAR04209', 'MAR06413', 'MAR06413', 'MAR06411', 'MAR06414', 'MAR06414'),
                         to=c('MAM01597m', 'MAM01587m', 'MAM01580m', 'MAM02183m', 'MAM02662m', 'MAM01306m', 'MAM01306m', 'MAM01306m', 'MAM02944m', 'MAM01597m', 'MAM02943m', 'MAM01597m', 'MAM02943m', 'MAM01255m', 'MAM02943m', 'MAM01862m', 'MAM01862m', 'MAM02439m', 'MAM02633m', 'MAM00765m', 'MAM02934m', 'MAM02984m', 'MAM02934m', 'MAM01701m', 'MAM02944m'),
                         arrows='to', dashes=FALSE, color='green', width=1)
tca_edges = rbind(tca_mrEdges, tca_rmEdges)

tca_fluxes_graph = function(rxn_fluxes, main=''){
  for(rxn in tca_rxns$id){
    tca_edges$width[tca_edges$from==rxn | tca_edges$to==rxn] = abs(rxn_fluxes[rxn])*5 + 1
    if(rxn_fluxes[rxn]==0){
      tca_edges$color[tca_edges$from==rxn | tca_edges$to==rxn] = 'grey'
    }
    else if(rxn%in%c('MAR04145', 'MAR04458', 'MAR04588', 'MAR04147', 'MAR04652', 'MAR08743', 'MAR04141') & rxn_fluxes[rxn]>0){
      tca_edges$color[tca_edges$from==rxn | tca_edges$to==rxn] = 'red'
    }
    else if(rxn%in%c('MAR04589', 'MAR04152', 'MAR04410', 'MAR06414') & rxn_fluxes[rxn]<0) tca_edges$color[tca_edges$from==rxn | tca_edges$to==rxn] = 'red'
  }
  visNetwork::visNetwork(rbind(tca_metabs, tca_rxns),
                         tca_edges, main=main)
}



# ------------------
# - Glutaminolysis -
# ------------------
glutaminolysis_rxns = data.frame(id=c('MAR03892', 'MAR03827', 'MAR04109', 'MAR03802', 'MAR03804', 'MAR09802', 'MAR03903', 'MAR03825', 'MAR03864', 'MAR05122'),
                                 label=c('MAR03892', 'MAR03827', 'MAR04109', 'MAR03802', 'MAR03804', 'MAR09802', 'MAR03903', 'MAR03825', 'MAR03864', 'MAR05122'),
                                 shape='box', shadow='FALSE', color.background=c(rep('lightgrey', 7), rep('thistle', 3)),
                                 color.border='black', font.color='black', font.size=20)
glutaminolysis_metabs = data.frame(id=c('MAM01975m', 'MAM01974m', 'MAM01306m', 'MAM01370m', 'MAM02633m', 'MAM01307m', 'MAM02819m', 'MAM01975c', 'MAM01974c', 'MAM01370c', 'MAM01369c', 'MAM00918m', 'MAM00918c'),
                                   label=c('glutamine', 'glutamate', 'AKG', 'aspartate', 'OAA', 'alanine', 'pyruvate', 'glutamine', 'glutamate', 'aspartate', 'asparagine', '3-sulfinoalanine', '3-sulfinoalanine'),
                                   shape='ellipse', shadow='TRUE', color.background=c(rep('gold', 7), rep('lightblue', 4), 'gold', 'lightblue'),
                                   color.border='black', font.color='black', font.size=20)
glutaminolysis_mrEdges = data.frame(from=c('MAM01975m', 'MAM01974m', 'MAM02633m', 'MAM01974m', 'MAM02819m', 'MAM01974m', 'MAM01974m', 'MAM01975c', 'MAM01370c', 'MAM01975c', 'MAM01370m', 'MAM01974c', 'MAM00918m', 'MAM01974c', 'MAM01974c'),
                                    to=c('MAR03892', 'MAR03827', 'MAR03827', 'MAR04109', 'MAR04109', 'MAR03802', 'MAR03804', 'MAR09802', 'MAR03903', 'MAR03903', 'MAR03825', 'MAR03825', 'MAR03864', 'MAR03864', 'MAR05122'),
                                    arrows='to', dashes=FALSE, color='green', width=1)
glutaminolysis_rmEdges = data.frame(from=c('MAR03892', 'MAR03827', 'MAR03827', 'MAR04109', 'MAR04109', 'MAR03802', 'MAR03804', 'MAR09802', 'MAR03903', 'MAR03903', 'MAR03825', 'MAR03825', 'MAR03864', 'MAR03864', 'MAR05122'),
                                    to=c('MAM01974m', 'MAM01306m', 'MAM01370m', 'MAM01307m', 'MAM01306m', 'MAM01306m', 'MAM01306m', 'MAM01974c', 'MAM01369c', 'MAM01974c', 'MAM01370c', 'MAM01974m', 'MAM00918c', 'MAM01974m', 'MAM01974m'),
                                    arrows='to', dashes=FALSE, color='green', width=1)
glutaminolysis_edges = rbind(glutaminolysis_mrEdges, glutaminolysis_rmEdges)

glutaminolysis_fluxes_graph = function(rxn_fluxes, main=''){
  for(rxn in glutaminolysis_rxns$id){
    glutaminolysis_edges$width[glutaminolysis_edges$from==rxn | glutaminolysis_edges$to==rxn] = abs(rxn_fluxes[rxn])*5 + 1
    if(rxn_fluxes[rxn]==0){
      glutaminolysis_edges$color[glutaminolysis_edges$from==rxn | glutaminolysis_edges$to==rxn] = 'grey'
    }
    else if(rxn%in%c('MAR03827', 'MAR04109', 'MAR03802', 'MAR03804') & rxn_fluxes[rxn]>0){
      glutaminolysis_edges$color[glutaminolysis_edges$from==rxn | glutaminolysis_edges$to==rxn] = 'red'
    }
  }
  visNetwork::visNetwork(rbind(glutaminolysis_metabs, glutaminolysis_rxns),
                         glutaminolysis_edges, main=main)
}



# ----------
# - OXPHOS -
# ----------
oxphos_rxns = data.frame(id=c('MAR06921', 'MAR06911', 'MAR06918', 'MAR06914', 'MAR13081', 'MAR06912', 'MAR06916'),
                         label=c('MAR06921', 'MAR06911', 'MAR06918', 'MAR06914', 'MAR13081', 'MAR06912', 'MAR06916'),
                         shape='box', shadow='FALSE', color.background='lightgrey',
                         color.border='black', font.color='black', font.size=20)
oxphos_metabs = data.frame(id=c('MAM02039m', 'MAM02553m', 'MAM02552m', 'MAM01803m', 'MAM01802m', 'MAM02630m', 'MAM02040m', 'MAM02631m', 'MAM02759m', 'MAM02751m', 'MAM02039i', 'MAM03103m', 'MAM03102m', 'MAM01824m', 'MAM01826m'),
                           label=c('H+', 'NADH', 'NAD+', 'FADH2', 'FAD', 'O2', 'H2O', 'O2-', 'PPi', 'Pi', 'H+', 'ubiquinone', 'ubiquinol', 'ferricytochrome C', 'ferrocytochrome C'),
                           shape=c(rep('circle', 11), rep('ellipse', 4)), shadow='TRUE', color.background=c(rep('gold',10), 'orange', rep('gold', 4)),
                           color.border='black', font.color='black', font.size=20)
oxphos_mrEdges = data.frame(from=c('MAM02039m', 'MAM02553m', 'MAM03103m', 'MAM01803m', 'MAM03103m', 'MAM02039m', 'MAM01824m', 'MAM03102m', 'MAM02039m', 'MAM02630m', 'MAM01826m', 'MAM02039m', 'MAM02630m', 'MAM01826m', 'MAM02040m', 'MAM02759m', 'MAM01285m', 'MAM02039i', 'MAM02751m'),
                            to=c('MAR06921', 'MAR06921', 'MAR06921', 'MAR06911', 'MAR06911', 'MAR06918', 'MAR06918', 'MAR06918', 'MAR06914', 'MAR06914', 'MAR06914', 'MAR13081', 'MAR13081', 'MAR13081', 'MAR06912', 'MAR06912', 'MAR06916', 'MAR06916', 'MAR06916'),
                            arrows='to', dashes=FALSE, color='green', width=1)
oxphos_rmEdges = data.frame(from=c('MAR06921', 'MAR06921', 'MAR06921', 'MAR06911', 'MAR06911', 'MAR06918', 'MAR06918', 'MAR06918', 'MAR06914', 'MAR06914', 'MAR06914', 'MAR13081', 'MAR13081', 'MAR13081', 'MAR13081', 'MAR06912', 'MAR06912', 'MAR06916', 'MAR06916', 'MAR06916'),
                            to=c('MAM02039i', 'MAM02552m', 'MAM03102m', 'MAM01802m', 'MAM03102m', 'MAM02039i', 'MAM01826m', 'MAM03103m', 'MAM02039i', 'MAM02040m', 'MAM01824m', 'MAM02039i', 'MAM02040m', 'MAM02631m', 'MAM01824m', 'MAM02039m', 'MAM02751m', 'MAM01371m', 'MAM02039m', 'MAM02040m'),
                            arrows='to', dashes=FALSE, color='green', width=1)
oxphos_edges = rbind(oxphos_mrEdges, oxphos_rmEdges)

oxphos_fluxes_graph = function(rxn_fluxes, main=''){
  for(rxn in oxphos_rxns$id){
    oxphos_edges$width[oxphos_edges$from==rxn | oxphos_edges$to==rxn] = abs(rxn_fluxes[rxn])*5 + 1
    if(rxn_fluxes[rxn]==0){
      oxphos_edges$color[oxphos_edges$from==rxn | oxphos_edges$to==rxn] = 'grey'
    }
    else if(rxn%in%c('MAR06911') & rxn_fluxes[rxn]<0) oxphos_edges$color[oxphos_edges$from==rxn | oxphos_edges$to==rxn] = 'red'
  }
  visNetwork::visNetwork(rbind(oxphos_metabs, oxphos_rxns),
                         oxphos_edges, main=main)
}




