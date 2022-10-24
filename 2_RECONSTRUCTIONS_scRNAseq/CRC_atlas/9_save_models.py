if __name__ == '__main__':
    from os.path import join, isdir
    from pandas import read_csv
    from os import listdir, mkdir
    from dill import load
    from gc import collect

    from cobra.io import write_sbml_model

    # Directories:
    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'
    CRCReconstruction_dir = join(base_dir, '2_RECONSTRUCTIONS_scRNAseq/CRC_atlas')
    CRCReconstructionNormalMatched_dir = join(CRCReconstruction_dir, 'NormalMatched')
    models_dir = join(base_dir, '0MODELS/CRC_atlas')

    # Files:
    pFBA_normal = join(CRCReconstructionNormalMatched_dir, 'FBA/Normal_Blood_pFBA_biomass.csv')
    CRCatlas_meta = join(CRCReconstructionNormalMatched_dir, 'metadata.csv')

    # Save Files:
    models_to_save = read_csv(CRCatlas_meta, index_col=0).index.tolist() # read_csv(pFBA_normal, index_col=0).columns.tolist()
    models_meta = read_csv(CRCatlas_meta, index_col=0)
    individuals = ['31', '32', '33', '35', 'KUL01', 'KUL19', 'KUL21', 'SMC01', 'SMC04', 'SMC06', 'SMC07',
                   'SMC08', 'SMC10']

    for indiv in individuals:
        print('\n', indiv)
        indiv_samples = [file for file in listdir(join(CRCReconstructionNormalMatched_dir, indiv))
                         if file.startswith('02_')]
        for samp in indiv_samples:
            samp_name = samp.replace('.obj', '').replace('02_', '')
            print('- ', samp_name)
            with open(join(CRCReconstructionNormalMatched_dir, indiv, samp), 'rb') as dump_file:
                temp_dump = load(dump_file)
                for cell_type, model in temp_dump.items():
                    model_name = '_'.join((indiv, samp_name, cell_type))
                    if model_name in models_to_save:
                        print('--', cell_type)
                        if not isdir(join(models_dir, cell_type)):
                            mkdir(join(models_dir, cell_type))
                        if not isdir(join(models_dir, cell_type, 'NormalMatched')):
                            mkdir(join(models_dir, cell_type, 'NormalMatched'))
                        if not isdir(join(models_dir, cell_type, 'Tumour')):
                            mkdir(join(models_dir, cell_type, 'Tumour'))
                        samp_state = models_meta.loc[indiv + '_' + samp_name + '_' + cell_type, 'state']
                        if samp_state == 'Tumour':
                            file_dir = join(models_dir, cell_type, 'Tumour')
                        else:
                            file_dir = join(models_dir, cell_type, 'NormalMatched')
                        file_name = file_dir + '/' + indiv + '_' + samp_name + '_' + cell_type + '.xml.gz'
                        write_sbml_model(model, file_name)
                        pass
                del temp_dump
                collect()
