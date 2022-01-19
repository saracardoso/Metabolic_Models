
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
