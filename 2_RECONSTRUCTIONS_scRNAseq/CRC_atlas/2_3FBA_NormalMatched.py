if __name__ == '__main__':
    from os.path import join
    from os import listdir
    from dill import load
    import json
    from pandas import read_csv, concat, DataFrame
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
    Evaluations upon medium changes
    '''

    # Run FBA in the different conditions:
    conditions = json.load(open(join(utilityData_dir, 'sampling_evals_Tcells.json')))
    CRCatlas_sampling = json.load(open(CRCatlas_sampling_file))
    predicted_fluxes = read_csv(join(CRCReconstructionNormalMatched_dir, '1_control_analysis/FBA/normal_FBA.csv'),
                                index_col=0)

    for condition_id in conditions.keys():
        print('\n', condition_id)
        res_fluxes = helper_fba_condition_normalmatched(conditions[condition_id], tcells_media_file,
                                                        'Blood_SMDB', CRCatlas_sampling,
                                                        CRCReconstructionNormalMatched_dir,
                                                        fba_results=predicted_fluxes)
        res_fluxes.to_csv(''.join((CRCReconstructionNormalMatched_dir, '/2_abnormal_medium/predicted_fluxes/',
                                   condition_id, '.csv')))

    '''
    Get reaction presence of models
    '''
    individuals = listdir(CRCReconstructionNormalMatched_dir)
    reactions_presence = DataFrame()
    CRCatlas_sampling = json.load(open(CRCatlas_sampling_file))
    for indiv in individuals:
        if search('[.]', indiv):
            next
        print('\n', indiv)
        # Get files that starts with 02_
        indiv_samples = [file for file in listdir(join(CRCReconstructionNormalMatched_dir, indiv)) if
                         file.startswith('02_')]
        for samp in indiv_samples:
            samp_name = samp.replace('.obj', '').replace('02_', '')
            print('- ', samp_name)
            with open(join(CRCReconstructionNormalMatched_dir, indiv, samp), 'rb') as dump_file:
                temp_dump = load(dump_file)
                for cell_type, model in temp_dump.items():
                    if cell_type in CRCatlas_sampling['NormalMatched'][indiv]['control'][samp_name] and cell_type in [
                        'Naive CD8 Tcells', 'Memory CD8 Tcells', 'Cytotoxic CD8 Tcells', 'Proliferative CD8 Tcells',
                        'Naive CD4 Tcells', 'Memory CD4 Tcells', 'Proliferative CD4 Tcells', 'IL17+ CD4 Tcells',
                        'Follicular CD4 Tcells', 'Regulatory CD4 Tcells']:
                        with_bounds = []
                        rxn_ids = []
                        for i in range(0, len(model.reactions)):
                            rxn_ids.append(model.reactions[i].id)
                            if model.reactions[i].bounds == (0, 0):
                                with_bounds.append(0)
                            else:
                                with_bounds.append(1)
                        reactions_presence['_'.join((indiv, samp_name, cell_type))] = with_bounds
                del temp_dump
                collect()
    reactions_presence = reactions_presence.set_axis(rxn_ids, axis='index')
    reactions_presence.to_csv(''.join((CRCReconstruction_dir, '/general/reaction_presence_Tcells.csv')))
