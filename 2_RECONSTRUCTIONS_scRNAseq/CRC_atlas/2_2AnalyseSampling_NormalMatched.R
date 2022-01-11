code_dir = './GENERAL/code/R'
for(file in list.files(code_dir, full.names=TRUE)) source(file)

# -----
# --- Create json file with names of models that 'passed' the sampling analysis [Control]
# --- [all samplings will be stored in this file]
# -----

CRCatlas_sampling = list()
CRCatlas_sampling$NormalMatched = list()

# Patient 31
CRCatlas_sampling$NormalMatched[['31']] = list()
CRCatlas_sampling$NormalMatched[['31']]$control = list()
CRCatlas_sampling$NormalMatched[['31']]$control$scrEXT001 = c('IL17+ CD4 Tcells', 'Memory CD4 Tcells',
                                                              'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['31']]$control$scrEXT002 = c('Cytotoxic CD8 Tcells',
                                                              'Memory CD4 Tcells',
                                                              'Proliferative CD8 Tcells',
                                                              'Follicular CD4 Tcells')
CRCatlas_sampling$NormalMatched[['31']]$control$scrEXT003 = c('Naive CD4 Tcells')

# Patient 32
CRCatlas_sampling$NormalMatched[['32']] = list()
CRCatlas_sampling$NormalMatched[['32']]$control = list()
CRCatlas_sampling$NormalMatched[['32']]$control$scrEXT009 = c('Cytotoxic CD8 Tcells',
                                                              'IL17+ CD4 Tcells', 'Memory CD8 Tcells',
                                                              'Naive CD8 Tcells',
                                                              'Prolifertive CD4 Tcells',
                                                              'Proliferative CD8 Tcells')
CRCatlas_sampling$NormalMatched[['32']]$control$scrEXT010 = c('Cytotoxic CD8 Tcells',
                                                              'Follicular CD4 Tcells',
                                                              'Memory CD4 Tcells',
                                                              'Memory CD8 Tcells', 'Naive CD4 Tcells',
                                                              'Naive CD8 Tcells',
                                                              'Prolifertive CD4 Tcells',
                                                              'Proliferative CD8 Tcells',
                                                              'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['32']]$control$scrEXT011 = c('Cytotoxic CD8 Tcells',
                                                              'Memory CD8 Tcells', 'Naive CD4 Tcells',
                                                              'Regulatory CD4 Tcells')

# Patient 33
CRCatlas_sampling$NormalMatched[['33']] = list()
CRCatlas_sampling$NormalMatched[['33']]$control = list()
CRCatlas_sampling$NormalMatched[['33']]$control$scrEXT012 = c('Memory CD8 Tcells',
                                                              'Prolifertive CD4 Tcells',
                                                              'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['33']]$control$scrEXT013 = c('Memory CD4 Tcells',
                                                              'Naive CD4 Tcells',
                                                              'Prolifertive CD4 Tcells')
CRCatlas_sampling$NormalMatched[['33']]$control$scrEXT014 = c('IL17+ CD4 Tcells',
                                                              'Memory CD4 Tcells',
                                                              'Memory CD8 Tcells', 'Naive CD4 Tcells')

# Patient 35
CRCatlas_sampling$NormalMatched[['35']] = list()
CRCatlas_sampling$NormalMatched[['35']]$control = list()
CRCatlas_sampling$NormalMatched[['35']]$control$scrEXT018 = c('Follicular CD4 Tcells',
                                                              'IL17+ CD4 Tcells', 'Memory CD4 Tcells',
                                                              'Memory CD8 Tcells',
                                                              'Naive CD4 Tcells',
                                                              'Prolifertive CD4 Tcells',
                                                              'Proliferative CD8 Tcells')
CRCatlas_sampling$NormalMatched[['35']]$control$scrEXT019 = c('Follicular CD4 Tcells',
                                                              'Memory CD4 Tcells',
                                                              'Memory CD8 Tcells',
                                                              'Naive CD4 Tcells',
                                                              'Reglatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['35']]$control$scrEXT020 = c('Memory CD8 Tcells', 'Naive CD4 Tcells',
                                                              'Regulatory CD4 Tcells')

# Patient KUL01
CRCatlas_sampling$NormalMatched[['KUL01']] = list()
CRCatlas_sampling$NormalMatched[['KUL01']]$control = list()
CRCatlas_sampling$NormalMatched[['KUL01']]$control[['KUL01-B']] = c('Cytotoxic CD8 Tcells',
                                                                 'Follicular CD4 Tcells',
                                                                 'Memory CD4 Tcells',
                                                                 'Memory CD8 Tcells',
                                                                 'Naive CD4 Tcells',
                                                                 'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['KUL01']]$control[['KUL01-N']] = c('Memory CD4 Tcells',
                                                                 'Memory CD8 Tcells',
                                                                 'Naive CD4 Tcells')
CRCatlas_sampling$NormalMatched[['KUL01']]$control[['KUL01-T']] = c('IL17+ CD4 Tcells',
                                                                 'Memory CD4 Tcells',
                                                                 'Memory CD8 Tcells',
                                                                 'Naive CD4 Tcells')

# Patient KUL19
CRCatlas_sampling$NormalMatched[['KUL19']] = list()
CRCatlas_sampling$NormalMatched[['KUL19']]$control = list()
CRCatlas_sampling$NormalMatched[['KUL19']]$control[['KUL19-B']] = c('Cytotoxic CD8 Tcells',
                                                                 'Follicular CD4 Tcells',
                                                                 'Memory CD4 Tcells',
                                                                 'Memory CD8 Tcells',
                                                                 'Naive CD4 Tcells',
                                                                 'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['KUL19']]$control[['KUL19-N']] = c('Memory CD4 Tcells',
                                                                 'Memory CD8 Tcells',
                                                                 'Naive CD4 Tcells')
CRCatlas_sampling$NormalMatched[['KUL19']]$control[['KUL19-T']] = c('IL17+ CD4 Tcells',
                                                                 'Memory CD4 Tcells',
                                                                 'Memory CD8 Tcells',
                                                                 'Naive CD4 Tcells')

# Patient KUL21
CRCatlas_sampling$NormalMatched[['KUL21']] = list()
CRCatlas_sampling$NormalMatched[['KUL21']]$control = list()
CRCatlas_sampling$NormalMatched[['KUL21']]$control[['KUL21-B']] = c('Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells')
CRCatlas_sampling$NormalMatched[['KUL21']]$control[['KUL21-N']] = c('IL17+ CD4 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['KUL21']]$control[['KUL21-T']] = c('Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')

# Patient SMC01
CRCatlas_sampling$NormalMatched[['SMC01']] = list()
CRCatlas_sampling$NormalMatched[['SMC01']]$control = list()
CRCatlas_sampling$NormalMatched[['SMC01']]$control[['SMC01-N']] = c('Cytotoxic CD8 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['SMC01']]$control[['SMC01-T']] = c('Cytotoxic CD8 Tcells',
                                                                    'Follicular CD4 Tcells',
                                                                    'IL17+ CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Proliferative CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')

# Patient SMC04
CRCatlas_sampling$NormalMatched[['SMC04']] = list()
CRCatlas_sampling$NormalMatched[['SMC04']]$control = list()
CRCatlas_sampling$NormalMatched[['SMC04']]$control[['SMC04-N']] = c('Cytotoxic CD8 Tcells',
                                                                    'IL17+ CD4 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Naive CD8 Tcells')
CRCatlas_sampling$NormalMatched[['SMC04']]$control[['SMC04-T']] = c('Cytotoxic CD8 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Naive CD8 Tcells',
                                                                    'Regulatory CD4 Tcells')

# Patient SMC06
CRCatlas_sampling$NormalMatched[['SMC06']] = list()
CRCatlas_sampling$NormalMatched[['SMC06']]$control = list()
CRCatlas_sampling$NormalMatched[['SMC06']]$control[['SMC06-N']] = c('Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Naive CD8 Tcells')
CRCatlas_sampling$NormalMatched[['SMC06']]$control[['SMC06-T']] = c('Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')

# Patient SMC07
CRCatlas_sampling$NormalMatched[['SMC07']] = list()
CRCatlas_sampling$NormalMatched[['SMC07']]$control = list()
CRCatlas_sampling$NormalMatched[['SMC07']]$control[['SMC07-N']] = c('Memory CD4 Tcells',
                                                                    'Naive CD4 Tcells')
CRCatlas_sampling$NormalMatched[['SMC07']]$control[['SMC07-T']] = c('Cytotoxic CD8 Tcells',
                                                                    'Follicular CD4 Tcells',
                                                                    'IL17+ CD4 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Naive CD8 Tcells',
                                                                    'Proliferative CD4 Tcells',
                                                                    'Proliferative CD8 Tcells',
                                                                    'Regulatory CD4 Tcells')

# Patient SMC08
CRCatlas_sampling$NormalMatched[['SMC08']] = list()
CRCatlas_sampling$NormalMatched[['SMC08']]$control = list()
CRCatlas_sampling$NormalMatched[['SMC08']]$control[['SMC08-N']] = c('Cytotoxic CD8 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['SMC08']]$control[['SMC08-T']] = c('Cytotoxic CD8 Tcells',
                                                                    'Follicular CD4 Tcells',
                                                                    'IL17+ CD4 Tcells',
                                                                    'Memory CD4 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Proliferative CD4 Tcells',
                                                                    'Regulatory CD4 Tcells')

# Patient SMC10
CRCatlas_sampling$NormalMatched[['SMC10']] = list()
CRCatlas_sampling$NormalMatched[['SMC10']]$control = list()
CRCatlas_sampling$NormalMatched[['SMC10']]$control[['SMC10-N']] = c('Cytotoxic CD8 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Regulatory CD4 Tcells')
CRCatlas_sampling$NormalMatched[['SMC10']]$control[['SMC10-T']] = c('Cytotoxic CD8 Tcells',
                                                                    'Follicular CD4 Tcells',
                                                                    'Memory CD8 Tcells',
                                                                    'Naive CD4 Tcells',
                                                                    'Proliferative CD4 Tcells',
                                                                    'Proliferative CD8 Tcells')

# Save the list into a json file:
jsonlite::write_json(CRCatlas_sampling, './2_RECONSTRUCTIONS_scRNAseq/CRC_atlas/CRCatlas_sampling.json')





