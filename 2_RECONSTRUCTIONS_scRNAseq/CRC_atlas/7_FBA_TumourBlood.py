if __name__ == '__main__':
    from os.path import join
    from os import listdir
    from dill import load
    import json
    from pandas import read_csv, concat, DataFrame
    from gc import collect
    from re import search

    '''
    Directories and files
    '''

    # Directories:
    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'
    utilityData_dir = join(base_dir, 'GENERAL/utility_data')
    CRCReconstruction_dir = join(base_dir, '2_RECONSTRUCTIONS_scRNAseq/CRC_atlas')
    CRCReconstructionNormalMatched_dir = join(CRCReconstruction_dir, 'NormalMatched')
    HumanGEM_dir = join(base_dir, '0MODELS/HumanGEM')

    # Files:
    tcells_media_file = join(utilityData_dir, 'media_Tcells_8percSerum.csv')
    CRCatlas_sampling_file = join(CRCReconstruction_dir, 'CRCatlas_sampling.json')
    HumanGem_file = join(HumanGEM_dir, 'HumanGEM-1.8.0_consistent.xml.gz')
    CRCatlas_meta = join(CRCReconstructionNormalMatched_dir, 'metadata.csv')

    '''
    Run FBA under normal medium conditions
    '''

    models_to_run = read_csv(CRCatlas_meta, index_col=0).index.tolist()
    individuals = ['31', '32', '33', '35', 'KUL01', 'KUL19', 'KUL21', 'SMC01', 'SMC04', 'SMC06', 'SMC07',
                   'SMC08', 'SMC10']
    predicted_fluxes = None
    models_names = []
    for indiv in individuals:
        print('\n', indiv)
        # Get files that starts with 02_
        indiv_samples = [file for file in listdir(join(CRCReconstructionNormalMatched_dir, indiv))
                         if file.startswith('02_')]
        for samp in indiv_samples:
            samp_name = samp.replace('.obj', '').replace('02_', '')
            print('- ', samp_name)
            with open(join(CRCReconstructionNormalMatched_dir, indiv, samp), 'rb') as dump_file:
                temp_dump = load(dump_file)
                for cell_type, model in temp_dump.items():
                    model_name = '_'.join((indiv, samp_name, cell_type))
                    if model_name in models_to_run:
                        print('--', cell_type)
                        # Get SMDB medium
                        media = read_csv(tcells_media_file, index_col='ID')
                        model.medium = media['Blood_tumour'].to_dict()
                        # Get objective
                        if cell_type in ['Proliferative CD4 Tcells', 'Proliferative CD8 Tcells']:
                            model.objective = {model.reactions.MAR13082: 1}
                        else:
                            model.objective = {model.reactions.MAR13082: 1, model.reactions.MAR06916: 1}
                        # Run pFBA
                        # model.solver = 'glpk_exact'
                        solFBA = model.optimize(objective_sense='maximize')
                        # solpFBA = pfba(model)
                        # Save predicted fluxes
                        models_names.append(model_name)
                        if predicted_fluxes is None:
                            predicted_fluxes = concat([solFBA.fluxes], axis=1)
                        else:
                            predicted_fluxes = concat([predicted_fluxes, solFBA.fluxes], axis=1)
                del temp_dump
                collect()
    predicted_fluxes.columns = models_names
    predicted_fluxes.to_csv(join(CRCReconstructionNormalMatched_dir, 'FBA/Tumour_Blood.csv'))
