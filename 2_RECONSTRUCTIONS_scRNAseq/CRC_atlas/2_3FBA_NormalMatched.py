if __name__ == '__main__':
    from os.path import join
    from os import listdir
    from dill import load
    import json
    from pandas import read_csv, concat
    from gc import collect
    from re import search

    from GENERAL.code.python.Validation import helper_fba_condition_normalmatched
    from cobra.flux_analysis import flux_variability_analysis

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

    '''
    Run FBA under normal medium conditions
    '''

    CRCatlas_sampling = json.load(open(CRCatlas_sampling_file))
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
                    if cell_type in CRCatlas_sampling['NormalMatched'][indiv]['control'][samp_name]:
                        print('--', cell_type)
                        # Get SMDB medium
                        media = read_csv(tcells_media_file, index_col='ID')
                        model.medium = media['Blood_SMDB'].to_dict()
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
                        models_names.append('_'.join((indiv, samp_name, cell_type)))
                        if predicted_fluxes is None:
                            predicted_fluxes = concat([solFBA.fluxes], axis=1)
                        else:
                            predicted_fluxes = concat([predicted_fluxes, solFBA.fluxes], axis=1)
                del temp_dump
                collect()
    predicted_fluxes.columns = models_names
    predicted_fluxes.to_csv(join(CRCReconstructionNormalMatched_dir, '1_control_analysis/FBA/normal_pFBA.csv'))

    '''
    Different condition evaluations
    '''

    # Run FVA:
    fva_rxns = ['MAR09517', 'MAR04171', 'MAR04193', 'MAR04210', 'MAR04358', 'MAR04421', 'MAR04573', 'MAR06627',
                'MAR09517']
    CRCatlas_sampling = json.load(open(CRCatlas_sampling_file))
    individuals = listdir(CRCReconstructionNormalMatched_dir)
    fva_fluxes = {}
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
                    if cell_type in CRCatlas_sampling['NormalMatched'][indiv]['control'][samp_name]:
                        print('--', cell_type)
                        # Get SMDB medium
                        media = read_csv(tcells_media_file, index_col='ID')
                        model.medium = media['Blood_SMDB'].to_dict()
                        # Get objective
                        if cell_type in ['Proliferative CD4 Tcells', 'Proliferative CD8 Tcells']:
                            model.objective = {model.reactions.MAR13082: 1}
                        else:
                            model.objective = {model.reactions.MAR13082: 1, model.reactions.MAR06916: 1}
                        # Run pFBA
                        # model.solver = 'glpk_exact'
                        solFVA = flux_variability_analysis(model, reaction_list=fva_rxns)
                        # Save predicted fluxes
                        fva_fluxes['_'.join((indiv, samp_name, cell_type))] = solFVA
                del temp_dump
                collect()
    with open(join(CRCReconstructionNormalMatched_dir, '1_control_analysis/FBA/FVA.json'), 'w') as f:
        json.dump(fva_fluxes, f)
    json.dump()


    # Run FBA in the different conditions:
    conditions = json.load(open(join(utilityData_dir, 'sampling_evals_Tcells.json')))
    CRCatlas_sampling = json.load(open(CRCatlas_sampling_file))
    predicted_fluxes = read_csv(join(CRCReconstructionNormalMatched_dir, '1_control_analysis/FBA/normal_FBA.csv'),
                                index_col=0)

    for condition_id in ['B10.1', 'B10.2', 'B10.3', 'B13', 'B14', 'B15.1', 'B15.2',
                         'B15.3', 'B18']: # 'B5', 'B7.1', 'B7.2', 'B7.3'
        # conditions.keys():
        print('\n', condition_id)
        res_fluxes = helper_fba_condition_normalmatched(conditions[condition_id], tcells_media_file,
                                                        'Blood_SMDB', CRCatlas_sampling,
                                                        CRCReconstructionNormalMatched_dir,
                                                        fba_results=predicted_fluxes)
        res_fluxes.to_csv(''.join((CRCReconstructionNormalMatched_dir, '/1_control_analysis/FBA/',
                                   condition_id, '.csv')))
