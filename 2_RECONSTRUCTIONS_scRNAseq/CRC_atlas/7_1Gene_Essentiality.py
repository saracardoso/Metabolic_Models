if __name__ == '__main__':
    from os.path import join
    from os import listdir
    from dill import load
    import json
    from pandas import read_csv, concat, Series
    from gc import collect
    from re import search

    from cobra.flux_analysis import single_gene_deletion
    from cobra import Solution
    from cobra.flux_analysis.room import room

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
    single_gene_deletion_file = join(utilityData_dir, 'single_gene_deletion.json')
    CRCatlas_meta = join(CRCReconstructionNormalMatched_dir, 'metadata.csv')

    '''
    Run single gene deletions
    '''

    models_to_run = read_csv(CRCatlas_meta, index_col=0).index.tolist()
    genes_rxns_toTest = json.load(open(single_gene_deletion_file))
    predicted_fluxes = read_csv(join(CRCReconstructionNormalMatched_dir, 'FBA/Normal_Blood.csv'),
                                index_col=0)
    individuals = ['31', '32', '33', '35', 'KUL01', 'KUL19', 'KUL21', 'SMC01', 'SMC04', 'SMC06', 'SMC07',
                   'SMC08', 'SMC10']
    gene_essentiality_growth = None
    gene_essentiality_status = None
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
                        model.medium = media['Blood_SMDB'].to_dict()
                        # Get objective
                        fba_orignal_fluxes = predicted_fluxes.loc[:, model_name]
                        if cell_type in ['Proliferative CD4 Tcells', 'Proliferative CD8 Tcells']:
                            model.objective = {model.reactions.MAR13082: 1}
                            fba_orig_value = fba_orignal_fluxes['MAR13082']
                        else:
                            model.objective = {model.reactions.MAR13082: 1, model.reactions.MAR06916: 1}
                            fba_orig_value = fba_orignal_fluxes['MAR13082'] + fba_orignal_fluxes['MAR06916']
                        # Get original FBA without knockouts:
                        fba_result = Solution(fba_orig_value, 'optimal', fluxes=fba_orignal_fluxes)
                        # Run gene essentiality
                        genes_toTest = [model.genes.get_by_id(gene) for gene in list(genes_rxns_toTest.keys())]
                        res_gess = single_gene_deletion(model, gene_list=genes_toTest,
                                                        method='linear room', solution=fba_result)
                        # Change frozenset indexes to strings:
                        new_indexes = []
                        for idx in res_gess.index:
                            x, *_ = idx
                            new_indexes.append(x)
                        res_gess = res_gess.set_axis(new_indexes, axis=0)
                        # Save results:
                        models_names.append(model_name)
                        if gene_essentiality_growth is None:
                            gene_essentiality_growth = concat([Series(res_gess['growth'].to_list())], axis=1)
                            gene_essentiality_status = concat([Series(res_gess['status'].to_list())], axis=1)
                        else:
                            gene_essentiality_growth = concat([gene_essentiality_growth,
                                                               Series(res_gess['growth'].to_list())], axis=1)
                            gene_essentiality_status = concat([gene_essentiality_status,
                                                               Series(res_gess['status'].to_list())], axis=1)
                del temp_dump
                collect()
    gene_essentiality_growth.columns = models_names
    gene_essentiality_growth = gene_essentiality_growth.set_axis(new_indexes, axis=0)
    gene_essentiality_status.columns = models_names
    gene_essentiality_status = gene_essentiality_status.set_axis(new_indexes, axis=0)

    gene_essentiality_growth.to_csv(''.join((CRCReconstructionNormalMatched_dir,
                                             '/Gene_essentiality/growth.csv')))
    gene_essentiality_status.to_csv(''.join((CRCReconstructionNormalMatched_dir,
                                             '/Gene_essentiality/status.csv')))

    '''
    Check glucose case:
    '''

    to_test = {}
    to_test['rxn'] = {}
    to_test['rxn']['SLC5A1'] = ['MAR01378', 'MAR06895', 'MAR08884']
    to_test['rxn']['GLUT1'] = ['MAR05029']
    to_test['rxn']['MAR05450'] = ['MAR05450']
    to_test['rxn']['all_uptakes'] = ['MAR01378', 'MAR06895', 'MAR08884', 'MAR05029', 'MAR05450']
    to_test['medium'] = {}
    to_test['medium']['glucose'] = ['MAR09034']
    to_test['medium']['glucose_sucrose'] = ['MAR09034', 'MAR09416']
    to_test['medium']['glucose_gcpool'] = ['MAR09034', 'MAR12043']
    to_test['medium']['glucose_sucrose_gcpool'] = ['MAR09034', 'MAR12043']
    to_test['invitro'] = {}
    to_test['invitro']['glucose'] = ['MAR09034']

    CRCatlas_sampling = json.load(open(CRCatlas_sampling_file))
    predicted_fluxes = read_csv(join(CRCReconstructionNormalMatched_dir, '1_control_analysis/FBA/normal_FBA.csv'),
                                index_col=0)
    individuals = listdir(CRCReconstructionNormalMatched_dir)
    for test_type in to_test.keys():
        for test_id in to_test[test_type].keys():
            print('\n', test_type)
            models_names = []
            predicted_fluxes2 = None
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
                                fba_orignal_fluxes = predicted_fluxes.loc[:, '_'.join((indiv, samp_name, cell_type))]
                                if cell_type in ['Proliferative CD4 Tcells', 'Proliferative CD8 Tcells']:
                                    model.objective = {model.reactions.MAR13082: 1}
                                    fba_orig_value = fba_orignal_fluxes['MAR13082']
                                else:
                                    model.objective = {model.reactions.MAR13082: 1, model.reactions.MAR06916: 1}
                                    fba_orig_value = fba_orignal_fluxes['MAR13082'] + fba_orignal_fluxes['MAR06916']
                                # Get original FBA without knockouts:
                                fba_result = Solution(fba_orig_value, 'optimal', fluxes=fba_orignal_fluxes)
                                # Get new predictions.
                                if test_type == 'rxn':
                                    with model as model_test:
                                        for rxn in to_test[test_type][test_id]:
                                            model_test.reactions.get_by_id(rxn).bounds = (0, 0)
                                        sol = room(model_test, solution=fba_result, linear=True)
                                else:
                                    with model as model_test:
                                        new_media = model_test.medium.copy()
                                        for metab in to_test[test_type][test_id]:
                                            new_media[metab] = 0
                                        model_test.medium = new_media.copy()
                                        sol = model_test.optimize(objective_sense='maximize')
                                # Save predicted fluxes
                                models_names.append('_'.join((indiv, samp_name, cell_type)))
                                if predicted_fluxes2 is None:
                                    predicted_fluxes2 = concat([sol.fluxes], axis=1)
                                else:
                                    predicted_fluxes2 = concat([predicted_fluxes2, sol.fluxes], axis=1)
                    del temp_dump
                    collect()
            predicted_fluxes2.columns = models_names
            predicted_fluxes2.to_csv(''.join((CRCReconstructionNormalMatched_dir,
                                              '/3_gene_essentiality/glucose_case/', test_id, '.csv')))
