
# ################################################# #
# INITIAL STEPS DONE TO OBTAIN BLOOD CONCENTRATIONS #
# ################################################# #

# Serum data was obtained from the following link:
# https://serummetabolome.ca/system/downloads/current/serum_metabolites.zip
# Release: 2019-01-17

serum_data_file = '/home/scardoso/Downloads/serum_metabolites.xml'
serum_data = XML::xmlParse(serum_data_file)
serum_data = XML::xmlToList(serum_data, simplify=T)

# Get normal blood concentrations into a dataframe (serum_data_normal_concentrations_all_info). others will store metabolites detected/
# expected but not quantified
serum_data_normal_concentrations_all_info = c()
others=c()
for(metabolite in serum_data){
  cat(metabolite$name, '\n')
  hmdb_id = metabolite$accession
  metab_name = metabolite$name
  pre = c(hmdb_id, metab_name)
  names(pre) = c('HMDB', 'Name')
  
  concentration_info=c()
  for(conc in metabolite$normal_concentrations){
    if(length(conc)>1){
      if(!is.null(conc$subject_condition)){
        if(conc$biospecimen == 'Blood' & !is.null(conc$concentration_value) & conc$subject_condition=='Normal' & !is.null(conc$subject_age)){
          x = as.data.frame(conc[1:4])
          if('comment'%in%names(conc)) x = cbind(x, conc$comment)
          else x = cbind(x, NA)
          colnames(x)[dim(x)[2]] = 'comment'
          x = c(pre, x)
          concentration_info = rbind(concentration_info, x)
        }
        else if(conc$biospecimen == 'Blood' & is.null(conc$concentration_value) & 'comment'%in%names(conc)){
          others = unique(rbind(others, c(hmdb_id, conc$comment)))
        }
      }
    }
  }
  serum_data_normal_concentrations_all_info = rbind(serum_data_normal_concentrations_all_info, concentration_info)
}
serum_data_normal_concentrations_all_info = as.data.frame(serum_data_normal_concentrations_all_info)

length(unique(serum_data_normal_concentrations_all_info$HMDB))

# Map HMDB IDs to HumanGEM IDs:
HumanGEM_mapping = as.data.frame(jsonlite::read_json('/home/scardoso/Downloads/humanGEMMetAssoc.JSON', simplifyVector = T))
HumanGEM_mapping_HG_HMDB = unique(HumanGEM_mapping[,c('mets', 'metHMDBID')])
HumanGEM_mapping_HG_HMDB = HumanGEM_mapping_HG_HMDB[grep('s$', HumanGEM_mapping_HG_HMDB$mets), ]
HumanGEM_mapping_HG_HMDB = HumanGEM_mapping_HG_HMDB[HumanGEM_mapping_HG_HMDB$metHMDBID!='', ]
for(i in 1:length(HumanGEM_mapping_HG_HMDB$metHMDBID)){
  hmdb = HumanGEM_mapping_HG_HMDB$metHMDBID[i]
  if(nchar(gsub('HMDB', '', hmdb))!=7){
    n_zeros = 7 - nchar(gsub('HMDB', '', hmdb))
    HumanGEM_mapping_HG_HMDB$metHMDBID[i] = paste('HMDB', paste(rep('0', n_zeros), collapse=''), gsub('HMDB', '', hmdb), sep='')
  }
}
HumanGEM_mapping_HG_HMDB[HumanGEM_mapping_HG_HMDB$mets=='pcholole_hs_s', 'metHMDBID'] = 'HMDB0002815'
HumanGEM_mapping_HG_HMDB[HumanGEM_mapping_HG_HMDB$mets=='pcholste_hs_s', 'metHMDBID'] = 'HMDB0010384'

HG_metab_ids = c()
for(hmdb in serum_data_normal_concentrations_all_info$HMDB){
  if(hmdb%in%HumanGEM_mapping_HG_HMDB$metHMDBID)
    HG_metab_ids = c(HG_metab_ids, HumanGEM_mapping_HG_HMDB[HumanGEM_mapping_HG_HMDB$metHMDBID==hmdb, 'mets'])
  else HG_metab_ids = c(HG_metab_ids, NA)
  if(length(HumanGEM_mapping_HG_HMDB[HumanGEM_mapping_HG_HMDB$metHMDBID==hmdb, 'mets'])>1)
    print(HumanGEM_mapping_HG_HMDB[HumanGEM_mapping_HG_HMDB$metHMDBID==hmdb, 'mets'])
}

colnames(HumanGEM_mapping_HG_HMDB)[2] = 'HMDB'

serum_data_normal_concentrations_all_info = plyr::join(HumanGEM_mapping_HG_HMDB, serum_data_normal_concentrations_all_info, by=c('HMDB'), type='full')


# Get data.frame with only HG ids with concentration data:
serum_concentrations = serum_data_normal_concentrations_all_info

# Parse the concentration_value variable into a numeric:
for(i in 1:dim(serum_concentrations)[1]){
  if(is.null(serum_concentrations$concentration_value[i][[1]]))
    serum_concentrations$concentration_value[i] = NA
  else{
    serum_concentrations$concentration_value[i] = strsplit(serum_concentrations$concentration_value[i][[1]], '[ ]*[(]')[[1]][1]
    serum_concentrations$concentration_value[i] = strsplit(serum_concentrations$concentration_value[i][[1]], '[ ]*[+]')[[1]][1]
    serum_concentrations$concentration_value[i] = strsplit(serum_concentrations$concentration_value[i][[1]], '[ ]*[±]')[[1]][1]
    if(length(grep('[-]', serum_concentrations$concentration_value[i][[1]]))!=0){
      splitted = strsplit(serum_concentrations$concentration_value[i][[1]], '[-]')[[1]]
      splitted = splitted[splitted!='']
      serum_concentrations$concentration_value[i] = mean(as.numeric(splitted))
    }
    else if(length(grep('^<', serum_concentrations$concentration_value[i][[1]]))!=0) serum_concentrations$concentration_value[i] = NA
    else if(length(grep('^>', serum_concentrations$concentration_value[i][[1]]))!=0) serum_concentrations$concentration_value[i] = NA
    else if(length(grep(',', serum_concentrations$concentration_value[i][[1]]))!=0) serum_concentrations$concentration_value[i] = NA
  }
  if(!is.na(serum_concentrations$concentration_value[i][[1]]) & is.na(as.numeric(serum_concentrations$concentration_value[i][[1]]))){
    cat(i, '\t')
    n_chars = nchar(serum_concentrations$concentration_value[i][[1]])
    serum_concentrations$concentration_value[i] = substr(serum_concentrations$concentration_value[i][[1]], 1, n_chars-1) 
  }
}
serum_concentrations$concentration_value = as.numeric(serum_concentrations$concentration_value)
serum_concentrations = serum_concentrations[-which(serum_concentrations$biospecimen=='Blood'&is.na(serum_concentrations$concentration_value)), 1:5]

# Average values from same metabolite:
concentrations_mets = c()
for(i in unique(na.omit(serum_concentrations$mets))){ # For mets not NA:
  concentration = mean(serum_concentrations$concentration_value[which(serum_concentrations$mets==i)])
  concentrations_mets = c(concentrations_mets, concentration)
}
names(concentrations_mets) = unique(na.omit(serum_concentrations$mets))
serum_concentrations_final = data.frame(HG=names(concentrations_mets), concentration=concentrations_mets)
colnames(HumanGEM_mapping_HG_HMDB)[1] = 'HG'
serum_concentrations_final = plyr::join(HumanGEM_mapping_HG_HMDB, serum_concentrations_final, by='HG', type = 'full')

concentrations_hmdbs = c()
for(i in unique(serum_concentrations[is.na(serum_concentrations$mets),'HMDB'])){ # For mets NA, go through HMDB ID:
  concentration = mean(serum_concentrations$concentration_value[which(serum_concentrations$HMDB==i)])
  concentrations_hmdbs = c(concentrations_hmdbs, concentration)
}
names(concentrations_hmdbs) = unique(serum_concentrations[is.na(serum_concentrations$mets),'HMDB'])
serum_concentrations_final = rbind(serum_concentrations_final,
                                   data.frame(HG=rep(NA, length(concentrations_hmdbs)),
                                              HMDB=names(concentrations_hmdbs),
                                              concentration=concentrations_hmdbs))
row.names(serum_concentrations_final) = as.character(1:dim(serum_concentrations_final)[1])

serum_concentrations_final = unique(plyr::join(serum_concentrations_final, serum_concentrations[,c('HMDB', 'Name')], by='HMDB', type='left'))
for(n in 1:dim(serum_concentrations_final)[1]) if(is.null(serum_concentrations_final$Name[n][[1]])) serum_concentrations_final$Name[n] = NA
serum_concentrations_final$Name = unlist(serum_concentrations_final$Name)

# Save serum_concentrations_final:
write.csv(serum_concentrations_final, '/home/scardoso/Downloads/serum_concentrations_final1.csv', row.names=F)

# This file was then manually curated (explanation of what was performed is elsewhere)



# ################################################################################## #
# OBTAIN FILE FOR GENE ID MAPPING OF CONSISTENT HUMANGEM MODEL, FOR SYMBOLS-ENSEMBLE #
# ################################################################################## #

alignment_mapping_file = './GENERAL/utility_data/genes_mapping_alignment.tsv'
humangem_mapping_file = './GENERAL/utility_data/genes_mapping_humangem.tsv'

alignment_mapping = read.table(alignment_mapping_file, header=TRUE, skip=1, sep='\t')
humangem_mapping = read.table(humangem_mapping_file, header=TRUE)

consistent_humangem_genes = read.table('./GENERAL/utility_data/HumanGEM-1.8.0_consistent_GENES.txt')$V1

alignment_mapping = alignment_mapping[,c(1,5)]
alignment_mapping = unique(alignment_mapping)
alignment_mapping = alignment_mapping[match(consistent_humangem_genes, alignment_mapping$Gene_stable_ID),]

humangem_mapping = humangem_mapping[,c(1,5)]
humangem_mapping = unique(humangem_mapping)
humangem_mapping = humangem_mapping[match(consistent_humangem_genes, humangem_mapping$genes),]

both_mappings = cbind(alignment_mapping, humangem_mapping$geneSymbols)
both_mappings[both_mappings[,2]!=both_mappings[,3],] # Ensemble genes with different gene symbols associated (same gene, but different names)

final_mapping = list()
for(gene in unique(c(both_mappings$Gene_name, both_mappings$`humangem_mapping$geneSymbols`))){
  if(gene%in%both_mappings$Gene_name) ensemble = both_mappings$Gene_stable_ID[both_mappings$Gene_name==gene]
  else ensemble = both_mappings$Gene_stable_ID[both_mappings$`humangem_mapping$geneSymbols`==gene]
  final_mapping[[gene]] = ensemble
}

jsonlite::write_json(final_mapping, './GENERAL/utility_data/genes_mapping.json')



# ##################################################################### #
# OBTAIN MAPPING OF REACTIONS IN SUBSYSTEMS AND SUBSYSTEMS OF REACTIONS #
# ##################################################################### #

HumanGEM_data_file = './GENERAL/utility_data/HumanGEM.json'

HumanGEM_data = jsonlite::read_json(HumanGEM_data_file, simplifyVector = T)$subsystems

# subsystems -> reactions
HumanGEM_rxns = HumanGEM_data$reactionList
names(HumanGEM_rxns) = HumanGEM_data$name
remove(HumanGEM_data)
invisible(gc())
HumanGEM_subsystems_rxns = HumanGEM_rxns
jsonlite::write_json(HumanGEM_subsystems_rxns,
                     './GENERAL/utility_data/subsystems_reactions_mapping.json')

# reactions -> subsystems
HumanGEM_rxns_subsystems = list()
for(subsystem in names(HumanGEM_subsystems_rxns)){
  rxns = HumanGEM_subsystems_rxns[[subsystem]]
  for(rxn in rxns){
    HumanGEM_rxns_subsystems[[rxn]] = c(HumanGEM_rxns_subsystems[[rxn]],
                                        subsystem)
  }
}
jsonlite::write_json(HumanGEM_rxns_subsystems,
                     './GENERAL/utility_data/reactions_subsystems_mapping.json')



# ############################################## #
# GATHER IMPORTANT REACTIONS FOR TCELLS ANALYSIS #
# ############################################## #

important_reactions_Tcells = list()

# Biomass:
important_reactions_Tcells$biomass = 'MAR13082'

# Reactions in subsystems:
fas = unique(c(HumanGEM_subsystems_rxns$`Fatty acid biosynthesis`,
             HumanGEM_subsystems_rxns$`Fatty acid biosynthesis (even-chain)`,
             HumanGEM_subsystems_rxns$`Fatty acid biosynthesis (odd-chain)`,
             HumanGEM_subsystems_rxns$`Fatty acid biosynthesis (unsaturated)`))
important_reactions_Tcells$subsystems[['Fatty acid biosynthesis']] = fas
important_reactions_Tcells$subsystems[['Fatty acid oxidation']] = HumanGEM_subsystems_rxns$`Fatty acid oxidation`
important_reactions_Tcells$subsystems[['Glycolysis / Gluconeogenesis']] = HumanGEM_subsystems_rxns$`Glycolysis / Gluconeogenesis`
important_reactions_Tcells$subsystems[['OXPHOS']] = HumanGEM_subsystems_rxns$`Oxidative phosphorylation`
important_reactions_Tcells$subsystems[['PPP']] = HumanGEM_subsystems_rxns$`Pentose phosphate pathway`
important_reactions_Tcells$subsystems[['TCA and glyoxylate/dicarboxylate metabolism']] = HumanGEM_subsystems_rxns$`Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism`
important_reactions_Tcells$subsystems[['Branched-chain amino acids (BCAA) catabolism']] = c('MAR03747', 'MAR03777', 'MAR06923', 'MAR03744', 'MAR03765', 'MAR03778')

# Glucose in PPP:
important_reactions_Tcells[['Glucose in PPP']] = c('MAR04306', 'MAR08971')

# Purines and pyrimidines precursor:
important_reactions_Tcells[['creation of purines precursor IMP']] = 'MAR04038'
important_reactions_Tcells[['creation of pyrimidines precursor orotate']] = 'MAR04577'

# Production of nucleotides:
important_reactions_Tcells[['nucleotide production']] = c('MAR07160')

# Protein synthesis:
important_reactions_Tcells[['protein synthesis']] = 'MAR13078'

# ATP Production:
important_reactions_Tcells[['protein synthesis']] = 'MAR06916'

# Uptakes:
important_reactions_Tcells$uptakes = list()
important_reactions_Tcells$uptakes[['GLUT1']] = c('MAR04858', 'MAR04973', 'MAR04996', 'MAR05027', 'MAR05029', 'MAR07675', 'MAR08846')
important_reactions_Tcells$uptakes[['SLC7A5']] = c('MAR05070', 'MAR05074', 'MAR05076', 'MAR05077', 'MAR05078', 'MAR05079',
                                                   'MAR05080', 'MAR05082', 'MAR05084', 'MAR05085', 'MAR05087', 'MAR05088',
                                                   'MAR05089', 'MAR05091', 'MAR05092', 'MAR05458', 'MAR05459', 'MAR05460',
                                                   'MAR05461', 'MAR05462', 'MAR05463', 'MAR05464', 'MAR05465', 'MAR05466',
                                                   'MAR05467', 'MAR05468', 'MAR05469', 'MAR05470', 'MAR05471', 'MAR05472',
                                                   'MAR05473', 'MAR05474', 'MAR05475', 'MAR05476', 'MAR05477', 'MAR05478',
                                                   'MAR05479', 'MAR05480', 'MAR05481', 'MAR05482', 'MAR05483', 'MAR05484',
                                                   'MAR05485', 'MAR05486', 'MAR05487', 'MAR05488', 'MAR05489', 'MAR05490',
                                                   'MAR05491', 'MAR05492', 'MAR05493', 'MAR05494', 'MAR05495', 'MAR05496',
                                                   'MAR05497', 'MAR05498', 'MAR05499', 'MAR05500', 'MAR05501', 'MAR05502',
                                                   'MAR05503', 'MAR05504', 'MAR05506', 'MAR05507', 'MAR05508', 'MAR05509',
                                                   'MAR05510', 'MAR05511', 'MAR05512', 'MAR05513', 'MAR05514', 'MAR05515',
                                                   'MAR05516', 'MAR05517', 'MAR05518', 'MAR05519', 'MAR05520', 'MAR05521',
                                                   'MAR05522', 'MAR05523', 'MAR05524', 'MAR05525', 'MAR05526', 'MAR05527',
                                                   'MAR05528', 'MAR05529', 'MAR05530', 'MAR05531', 'MAR05532', 'MAR05533',
                                                   'MAR05534', 'MAR05535', 'MAR05536', 'MAR05537', 'MAR05538', 'MAR05539',
                                                   'MAR05540', 'MAR05541', 'MAR05542', 'MAR05543', 'MAR05544', 'MAR05545',
                                                   'MAR05546', 'MAR05547', 'MAR05548', 'MAR05549', 'MAR05551', 'MAR05552',
                                                   'MAR05553', 'MAR05554', 'MAR05555', 'MAR05556', 'MAR05557', 'MAR05558',
                                                   'MAR05559', 'MAR05560', 'MAR05561', 'MAR05562', 'MAR05563', 'MAR05564',
                                                   'MAR05565', 'MAR05566', 'MAR05567', 'MAR05568', 'MAR05569', 'MAR05570',
                                                   'MAR05571', 'MAR05572', 'MAR05573', 'MAR05574', 'MAR05575', 'MAR05576',
                                                   'MAR05577', 'MAR05578', 'MAR05579', 'MAR05580', 'MAR05581', 'MAR05582',
                                                   'MAR05583', 'MAR05584', 'MAR05585', 'MAR05586', 'MAR05587', 'MAR05588',
                                                   'MAR05589', 'MAR05590', 'MAR05591', 'MAR05592', 'MAR05593', 'MAR05594',
                                                   'MAR05595', 'MAR05596', 'MAR05597', 'MAR05598', 'MAR05599', 'MAR05601',
                                                   'MAR05602', 'MAR05605', 'MAR05606', 'MAR05607', 'MAR05608', 'MAR05610',
                                                   'MAR05611', 'MAR05612', 'MAR05613', 'MAR05614', 'MAR05615', 'MAR05616',
                                                   'MAR05620', 'MAR05621', 'MAR05622', 'MAR05623', 'MAR05625', 'MAR05626',
                                                   'MAR05628', 'MAR05630', 'MAR05631', 'MAR05632', 'MAR05635', 'MAR05636',
                                                   'MAR05637', 'MAR05638', 'MAR05640', 'MAR05641', 'MAR05643', 'MAR05644',
                                                   'MAR05645', 'MAR05646', 'MAR05647', 'MAR05648', 'MAR05649', 'MAR05650',
                                                   'MAR05651', 'MAR05652', 'MAR05653', 'MAR05654', 'MAR05655', 'MAR05656',
                                                   'MAR05657', 'MAR05658', 'MAR05659', 'MAR05660', 'MAR05661', 'MAR05662',
                                                   'MAR05663', 'MAR05664', 'MAR05665', 'MAR05666', 'MAR05667', 'MAR05668',
                                                   'MAR05669', 'MAR05670', 'MAR05671', 'MAR05672', 'MAR05673', 'MAR05674',
                                                   'MAR05675', 'MAR05676', 'MAR06859', 'MAR09736')
important_reactions_Tcells$uptakes[['SLC1A5']] = c('MAR05305', 'MAR05308', 'MAR05603', 'MAR05604', 'MAR05609', 'MAR05612',
                                                   'MAR05614', 'MAR05618', 'MAR05619', 'MAR05624', 'MAR05627', 'MAR05629',
                                                   'MAR05633', 'MAR05634', 'MAR05639', 'MAR05642', 'MAR05644', 'MAR05708',
                                                   'MAR05709', 'MAR05710', 'MAR05717', 'MAR05719', 'MAR05753', 'MAR05754',
                                                   'MAR05755', 'MAR05760', 'MAR05764', 'MAR05783', 'MAR05784', 'MAR05785',
                                                   'MAR05790', 'MAR05793', 'MAR06377', 'MAR11808')
important_reactions_Tcells$uptakes[['SLC3A2']] = c('MAR04931', 'MAR05070', 'MAR05073', 'MAR05074', 'MAR05076', 'MAR05077',
                                                   'MAR05078', 'MAR05079', 'MAR05080', 'MAR05082', 'MAR05084', 'MAR05085',
                                                   'MAR05088', 'MAR05089', 'MAR05091', 'MAR05092', 'MAR05315', 'MAR05458',
                                                   'MAR05460', 'MAR05465', 'MAR05470', 'MAR05473', 'MAR05487', 'MAR05492',
                                                   'MAR05497', 'MAR05516', 'MAR05578', 'MAR05873', 'MAR07643', 'MAR07644',
                                                   'MAR08733', 'MAR08734', 'MAR09190', 'MAR09191', 'MAR00072', 'MAR06859',
                                                   'MAR06863', 'MAR06897', 'MAR07102', 'MAR07105', 'MAR07151', 'MAR07153',
                                                   'MAR07189', 'MAR11804', 'MAR11810')
important_reactions_Tcells$uptakes[['SLC38A2']] = c('MAR05305', 'MAR05308', 'MAR05311', 'MAR05315', 'MAR05319', 'MAR05320',
                                                    'MAR05322', 'MAR05324', 'MAR05326', 'MAR05333', 'MAR06375', 'MAR06380',
                                                    'MAR06381', 'MAR06382', 'MAR06384', 'MAR11791')
important_reactions_Tcells$uptakes[['CD36']] = c('MAR00030', 'MAR00155', 'MAR00164', 'MAR00167', 'MAR00170', 'MAR00173',
                                                 'MAR00176', 'MAR00179', 'MAR00183', 'MAR00187', 'MAR00191', 'MAR00195',
                                                 'MAR00199', 'MAR00203', 'MAR00208', 'MAR00212', 'MAR00216', 'MAR00225',
                                                 'MAR00232', 'MAR00236', 'MAR00240', 'MAR00244', 'MAR00248', 'MAR00254',
                                                 'MAR00258', 'MAR00262', 'MAR00266', 'MAR00270', 'MAR00274', 'MAR00278',
                                                 'MAR00282', 'MAR00288', 'MAR00292', 'MAR00296', 'MAR00300', 'MAR00304',
                                                 'MAR00308', 'MAR00312', 'MAR00318', 'MAR00322', 'MAR00326', 'MAR00330',
                                                 'MAR00336', 'MAR00340', 'MAR00344', 'MAR00348', 'MAR00352', 'MAR00356',
                                                 'MAR00360', 'MAR00364', 'MAR00368', 'MAR00372', 'MAR00376', 'MAR00380',
                                                 'MAR00384', 'MAR00388', 'MAR00392', 'MAR00396', 'MAR00400', 'MAR00404',
                                                 'MAR00408', 'MAR00412', 'MAR00416', 'MAR00420', 'MAR00424', 'MAR00428',
                                                 'MAR00432', 'MAR00436', 'MAR11815')
important_reactions_Tcells$uptakes[['fatty acids']] = c('MAR00240', 'MAR00428', 'MAR00380', 'MAR00288', 'MAR00254', 'MAR00432',
                                                        'MAR00424', 'MAR00274', 'MAR00420', 'MAR00372', 'MAR00270', 'MAR00203',
                                                        'MAR00199', 'MAR00416', 'MAR00368', 'MAR00388', 'MAR00436', 'MAR00392',
                                                        'MAR00384', 'MAR00232', 'MAR00300', 'MAR00296', 'MAR00244', 'MAR00376',
                                                        'MAR00364', 'MAR00360', 'MAR00412', 'MAR00408', 'MAR00312', 'MAR00340',
                                                        'MAR00322', 'MAR00318', 'MAR00292', 'MAR00258', 'MAR00404', 'MAR00282',
                                                        'MAR00266', 'MAR00400', 'MAR00308', 'MAR00187', 'MAR00330', 'MAR00396',
                                                        'MAR00348', 'MAR00236', 'MAR00304', 'MAR00195', 'MAR00336', 'MAR00278',
                                                        'MAR00262', 'MAR00356', 'MAR00216', 'MAR00225', 'MAR00212', 'MAR00208',
                                                        'MAR00248', 'MAR00352', 'MAR00326', 'MAR00191', 'MAR00344')
important_reactions_Tcells$uptakes[['ASCT2/SLC1A5']] = c('MAR05604', 'MAR05609', 'MAR05633', 'MAR05708', 'MAR05639', 'MAR05710')
important_reactions_Tcells$uptakes[['proline']] = c('MAR05567', 'MAR05566', 'MAR05565', 'MAR05564', 'MAR05563', 'MAR05557',
                                                    'MAR05551', 'MAR05544', 'MAR05536', 'MAR05527', 'MAR05517', 'MAR05506',
                                                    'MAR05494', 'MAR05481', 'MAR05467', 'MAR05737', 'MAR05738', 'MAR05739',
                                                    'MAR05740', 'MAR05741', 'MAR05742', 'MAR05743', 'MAR05744', 'MAR05745',
                                                    'MAR05746', 'MAR05747', 'MAR05748', 'MAR05749', 'MAR05750', 'MAR05751',
                                                    'MAR05877', 'MAR05878', 'MAR05879', 'MAR05880', 'MAR05881', 'MAR05322',
                                                    'MAR06373', 'MAR07194', 'MAR08938', 'MAR09617', 'MAR08988', 'MAR08929')
important_reactions_Tcells$uptakes[['arginine']] = c('MAR05068', 'MAR05923', 'MAR05924', 'MAR05925', 'MAR05926', 'MAR05927',
                                                     'MAR05928', 'MAR05929', 'MAR05930', 'MAR05931', 'MAR05932', 'MAR05933',
                                                     'MAR05934', 'MAR05935', 'MAR05936', 'MAR05937', 'MAR05938', 'MAR08629',
                                                     'MAR09592', 'MAR09871', 'MAR11804', 'MAR11810', 'MAR05578', 'MAR05582',
                                                     'MAR05583', 'MAR05584', 'MAR09182')
important_reactions_Tcells$uptakes[['CD98']] = c('MAR00072', 'MAR04931', 'MAR05070', 'MAR05073', 'MAR05074', 'MAR05076',
                                                 'MAR05077', 'MAR05078', 'MAR05079', 'MAR05080', 'MAR05082', 'MAR05084',
                                                 'MAR05085', 'MAR05088', 'MAR05089', 'MAR05091', 'MAR05092', 'MAR05315',
                                                 'MAR05458', 'MAR05460', 'MAR05465', 'MAR05470', 'MAR05473', 'MAR05487',
                                                 'MAR05492', 'MAR05497', 'MAR05516', 'MAR05578', 'MAR05873', 'MAR06859',
                                                 'MAR06863', 'MAR06897', 'MAR07102', 'MAR07105', 'MAR07151', 'MAR07153',
                                                 'MAR07189', 'MAR07643', 'MAR07644', 'MAR08733', 'MAR08734', 'MAR09190',
                                                 'MAR09191', 'MAR11804', 'MAR11810')
important_reactions_Tcells$uptakes[['glutamine']] = c('MAR05308', 'MAR05310', 'MAR05617', 'MAR05618', 'MAR05619', 'MAR05620',
                                                      'MAR05621', 'MAR05622', 'MAR05623', 'MAR05624', 'MAR05625', 'MAR05626',
                                                      'MAR05627', 'MAR05628', 'MAR05629', 'MAR05630', 'MAR05631', 'MAR05837',
                                                      'MAR05838', 'MAR05839', 'MAR05840', 'MAR05841', 'MAR06380', 'MAR06896',
                                                      'MAR06897', 'MAR09887', 'MAR05078', 'MAR05459', 'MAR05473', 'MAR05487',
                                                      'MAR05488', 'MAR05489', 'MAR05490', 'MAR05491', 'MAR05492', 'MAR05493',
                                                      'MAR05494', 'MAR05495', 'MAR05496', 'MAR05497', 'MAR05498', 'MAR05499',
                                                      'MAR07643', 'MAR08733')
important_reactions_Tcells$uptakes[['acetate']] = c('MAR04934', 'MAR04935', 'MAR06001', 'MAR06008', 'MAR06015', 'MAR06022',
                                                    'MAR06031', 'MAR06032', 'MAR06033', 'MAR06034', 'MAR06035', 'MAR06036')

# Other transport reactions:
important_reactions_Tcells$other_tranports[['pyruvate_into_mitochondria']] = c('MAR04926', 'MAR05426')

# Other enzymes:
important_reactions_Tcells[['HK2']] = c('MAR04319', 'MAR04394', 'MAR04486', 'MAR04487', 'MAR04488', 'MAR04489', 'MAR04490',
                                        'MAR04492', 'MAR04493', 'MAR04494', 'MAR04495', 'MAR04496', 'MAR07746')
important_reactions_Tcells[['PFK2']] = c('MAR04379')
important_reactions_Tcells[['MCT1']] = c('MAR01517', 'MAR02151')
important_reactions_Tcells[['LDHA']] = c('MAR01453', 'MAR04280', 'MAR04281', 'MAR04287', 'MAR04388', 'MAR11476', 'MAR11479')
important_reactions_Tcells[['CS']] = c('MAR04145')
important_reactions_Tcells[['OGDH']] = c('MAR04209', 'MAR04239', 'MAR04599', 'MAR05297', 'MAR06411', 'MAR06413')
important_reactions_Tcells[['CytC']] = c('MAR06918')
important_reactions_Tcells[['ATP5A']] = c('MAR06916')
important_reactions_Tcells[['ACC']] = c('MAR00030', 'MAR04156', 'MAR07672', 'MAR07673', 'MAR04295')
important_reactions_Tcells[['GLS']] = c('MAR03892')
important_reactions_Tcells[['GOT2']] = c('MAR01616', 'MAR03827', 'MAR03865', 'MAR03912', 'MAR04786', 'MAR06725', 'MAR06768',
                                         'MAR06780', 'MAR08066', 'MAR08067', 'MAR08098', 'MAR08684')
important_reactions_Tcells[['CPT1A']] = c('MAR00027', 'MAR00030', 'MAR00103', 'MAR00108', 'MAR00159', 'MAR00452', 'MAR01488',
                                          'MAR02121', 'MAR02122', 'MAR02123', 'MAR02591', 'MAR02594', 'MAR02599', 'MAR02602',
                                          'MAR02605', 'MAR02608', 'MAR02611', 'MAR02614', 'MAR02620', 'MAR02623', 'MAR02626',
                                          'MAR02633', 'MAR02636', 'MAR02642', 'MAR02648', 'MAR02651', 'MAR02654', 'MAR02660',
                                          'MAR02663', 'MAR02666', 'MAR02669', 'MAR02672', 'MAR02675', 'MAR02678', 'MAR02681',
                                          'MAR02684', 'MAR02687', 'MAR02690', 'MAR02693', 'MAR02699', 'MAR02702', 'MAR02705',
                                          'MAR02708', 'MAR02711', 'MAR02714', 'MAR02717', 'MAR02720', 'MAR02726', 'MAR02730',
                                          'MAR02733', 'MAR02736', 'MAR02739', 'MAR02742', 'MAR02746', 'MAR02752', 'MAR02755',
                                          'MAR02759', 'MAR02762', 'MAR02768', 'MAR02771', 'MAR02774', 'MAR02777', 'MAR02906',
                                          'MAR02909', 'MAR02912', 'MAR02915', 'MAR02918', 'MAR02921', 'MAR02924', 'MAR02927',
                                          'MAR02930', 'MAR03028', 'MAR03032', 'MAR04457', 'MAR04462', 'MAR04475', 'MAR04491',
                                          'MAR04504', 'MAR04506', 'MAR04508', 'MAR04522', 'MAR04562', 'MAR04569', 'MAR04571',
                                          'MAR04581', 'MAR04598', 'MAR04610', 'MAR04613', 'MAR04620', 'MAR04624', 'MAR04634',
                                          'MAR04636', 'MAR04638', 'MAR04639', 'MAR04671', 'MAR04707', 'MAR04711', 'MAR04719',
                                          'MAR04745', 'MAR04749', 'MAR04759', 'MAR05104', 'MAR05434', 'MAR05456', 'MAR05991',
                                          'MAR06278', 'MAR06598', 'MAR06696')
important_reactions_Tcells[['PCX']] = c('MAR04143')
important_reactions_Tcells[['PDH']] = c('MAR04137')

# Specific reactions:
important_reactions_Tcells$specific[['pyruvate_to_acetylCoA']] = c('MAR04137')
important_reactions_Tcells$specific[['PEP_to_pyruvate']] = c('MAR04358')
important_reactions_Tcells$specific[['pyruvate_to_lactate']] = c('MAR04388')
important_reactions_Tcells$specific[['3-PG consumption in serie biosynthesis']] = c('MAR03839')

# Others:
important_reactions_Tcells[['ornithine_production']] = c('MAR08426', 'MAR03807', 'MAR03809', 'MAR03816', 'MAR04212',
                                                         'MAR08425', 'MAR04582')
important_reactions_Tcells[['putrescine_production']] = c('MAR06961', 'MAR04212', 'MAR04424')

# Save information in a json file:
jsonlite::write_json(important_reactions_Tcells,
                     './GENERAL/utility_data/important_reactions_Tcells.json')

