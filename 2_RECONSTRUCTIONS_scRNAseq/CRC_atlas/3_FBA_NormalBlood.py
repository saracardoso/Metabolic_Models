if __name__ == '__main__':
    from os.path import join
    from os import listdir
    from dill import load
    import json
    from pandas import read_csv, concat, DataFrame
    from gc import collect
    from re import search

    from cobra.flux_analysis import pfba

    from GENERAL.code.python.Validation import helper_fba_condition_normalmatched

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
    CRCatlas_meta = join(CRCReconstructionNormalMatched_dir, 'metadata.csv')
    HumanGem_file = join(HumanGEM_dir, 'HumanGEM-1.8.0_consistent.xml.gz')
    pFBA_normal = join(CRCReconstructionNormalMatched_dir, 'FBA/Normal_Blood_pFBA_biomass.csv')

    '''
    Run FBA under normal medium conditions
    '''

    models_to_run = read_csv(CRCatlas_meta, index_col=0).index.tolist()
    individuals = listdir(CRCReconstructionNormalMatched_dir)
    predicted_fluxes = None
    models_names = []
    for indiv in individuals:
        if search('[.]', indiv):
            next
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
                        model.medium = media['Blood_SMDB'].to_dict()
                        # Get objective
                        # -- This next section is removed when all models are to be set to only maximise biomass
                        if cell_type in ['Proliferative CD4 Tcells', 'Proliferative CD8 Tcells']:
                            model.objective = {model.reactions.MAR13082: 1}
                        else:
                            model.objective = {model.reactions.MAR13082: 1, model.reactions.MAR06916: 1}
                        #model.objective = {model.reactions.MAR13082: 1}
                        # --
                        # Run pFBA
                        # model.solver = 'glpk_exact'
                        #solFBA = model.optimize(objective_sense='maximize')
                        try:
                            solFBA = pfba(model)
                        except:
                            print('Infeasible')
                        else:
                            # Save predicted fluxes
                            models_names.append(model_name)
                            if predicted_fluxes is None:
                                predicted_fluxes = concat([solFBA.fluxes], axis=1)
                            else:
                                predicted_fluxes = concat([predicted_fluxes, solFBA.fluxes], axis=1)
                del temp_dump
                collect()
    predicted_fluxes.columns = models_names
    # predicted_fluxes.to_csv(join(CRCReconstructionNormalMatched_dir, 'FBA/Normal_Blood.csv'))
    predicted_fluxes.to_csv(join(CRCReconstructionNormalMatched_dir, 'FBA/Normal_Blood_pFBA.csv'))

    '''
    Evaluations upon medium changes
    '''

    # Run FBA in the different conditions:
    conditions = json.load(open(join(utilityData_dir, 'evals_Tcells.json')))
    #models_to_run = read_csv(CRCatlas_meta, index_col=0).index.tolist()
    models_to_run = read_csv(pFBA_normal, index_col=0).columns.tolist()
    predicted_fluxes = read_csv(join(CRCReconstructionNormalMatched_dir, 'FBA/Normal_Blood.csv'),
                                index_col=0)
    individuals = ['31', '32', '33', '35', 'KUL01', 'KUL19', 'KUL21', 'SMC01', 'SMC04', 'SMC06', 'SMC07',
                   'SMC08', 'SMC10']

    for condition_id in ['B2', 'B4', 'B8.2_and_B9', 'B17', 'B8.5']: #['B6.1', 'B6.2', 'B8.1', 'B8.2_and_B9', 'B8.4', 'B8.5', 'B11', 'B12.1',
                         #'B12.2', 'B12.3', 'B16', 'B17']:#['B1', 'B2', 'B4', 'B6.1', 'B6.2', 'B8.1', 'B8.2_and_B9', 'B8.4', 'B8.5', 'B11', 'B12.1',
                         #'B12.2', 'B12.3', 'B16', 'B17']:
        print('\n', condition_id)
        res_fluxes = helper_fba_condition_normalmatched(conditions[condition_id], tcells_media_file,
                                                        'Blood_SMDB', models_to_run,
                                                        CRCReconstructionNormalMatched_dir,
                                                        individuals)#,
                                                        #fba_results=predicted_fluxes)
        res_fluxes.to_csv(''.join((CRCReconstructionNormalMatched_dir, '/FBA/medium_changes/',
                                   condition_id, '_biomass.csv')))

    # Chloride case:
    chloride_condition = {'change_media_bounds': {'MAR09150': 0}, 'change_internal_bounds': None}
    x = helper_fba_condition_normalmatched(chloride_condition, tcells_media_file, 'Blood_SMDB', models_to_run,
                                           CRCReconstructionNormalMatched_dir, individuals)
    x.to_csv(''.join((CRCReconstructionNormalMatched_dir, '/FBA/medium_changes/chloride_biomass.csv')))

    # Chloride case:
    sodium_condition = {'change_media_bounds': {'MAR09150': 0}, 'change_internal_bounds': None}
    x = helper_fba_condition_normalmatched(chloride_condition, tcells_media_file, 'Blood_SMDB', models_to_run,
                                           CRCReconstructionNormalMatched_dir, individuals)
    x.to_csv(''.join((CRCReconstructionNormalMatched_dir, '/FBA/medium_changes/chloride_biomass.csv')))
