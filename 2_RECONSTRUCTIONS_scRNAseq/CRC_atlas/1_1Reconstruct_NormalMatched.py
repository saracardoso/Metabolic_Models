if __name__ == '__main__':
    import numpy as np
    import pandas as pd
    from os.path import join, isdir
    from os import listdir, mkdir
    from pandas import read_csv
    from dill import dump, load

    from cobra.io import read_sbml_model

    from GENERAL.code.python.ReconstructModel import ReconstructModel
    from GENERAL.code.python.EvaluateModel import GapFillModel

    all_cell_types = ['Normal Epithelial cells', 'Tumour Epithelial cells', 'CAFs', 'Fibroblasts', 'VSMCs',
                      'Vascular ECs', 'Lymphatic ECs', 'Pericytes', 'Myofibroblasts', 'Enteric glia cells',
                      'Anti-inflammatory macro/mono', 'Pro-inflammatory macro/mono', 'Mast cells', 'cDCs', 'pDCs',
                      'Naive CD4 Tcells', 'Memory CD4 Tcells', 'Proliferative CD4 Tcells', 'Regulatory CD4 Tcells',
                      'IL17+ CD4 Tcells', 'IL22+ CD4 Tcells', 'Follicular CD4 Tcells', 'Naive CD8 Tcells',
                      'Cytotoxic CD8 Tcells', 'Memory CD8 Tcells', 'Proliferative CD8 Tcells', 'CXCL13+ CD8 Tcells',
                      'gdTcells', 'Double-Negative Tcells', 'CD8aa IELs', 'LTi cells', 'NK cells', 'NKT cells',
                      'Naive Bcells', 'Memory Bcells', 'Proliferative Bcells', 'Plasma cells']

    '''
    Directories and files
    '''

    # Directories:
    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'
    HumanGEM_dir = join(base_dir, '0MODELS/HumanGEM')
    utilityData_dir = join(base_dir, 'GENERAL/utility_data')
    CRCatlasNormalMatched_dir = join(base_dir, '0Data/scRNAseq/CRC_atlas/expression_data/CRC/normal_matched')
    CRCatlasReconstruction_dir = join(base_dir, '2_RECONSTRUCTIONS_scRNAseq/CRC_atlas')

    # Files:
    HumanGem_file = join(HumanGEM_dir, 'HumanGEM-1.8.0_consistent.xml.gz')
    mediaOpenBounds_file = join(utilityData_dir, 'media_openBounds.csv')
    genesMapping_file = join(utilityData_dir, 'genes_mapping.json')

    '''
    Read HumanGEM model
    '''

    print('\nReading HumanGEM-1.8.0...')
    HumanGEM = read_sbml_model(HumanGem_file)


    '''
    Get reactions present in media
    '''

    print('Getting media reactions...')
    media_rxns = read_csv(mediaOpenBounds_file, index_col='ID').index.tolist()

    '''
    Reconstruct model
    '''

    print('\nReconstructing models...')
    individuals = listdir(CRCatlasNormalMatched_dir)
    reconstructions = {}
    for individual in ['SMC07', 'SMC08', 'SMC10']:#individuals:
        print('.. ', individual)
        reconstructions[individual] = {}
        samples = listdir(join(CRCatlasNormalMatched_dir, individual))
        samples = [samp.replace('.csv', '') for samp in samples]
        nReactions_reconstruction = np.empty((len(all_cell_types), len(samples)))
        nReactions_reconstruction[:] = np.NaN
        nReactions_reconstruction = pd.DataFrame(nReactions_reconstruction, index=all_cell_types, columns=samples)
        for sample in samples:
            samp_fullpath = join(CRCatlasNormalMatched_dir, individual, ''.join((sample, '.csv')))
            print('... ', sample)
            # Reconstruct ...
            print('\n.... Reconstruction')
            reconstructions[individual][sample] = ReconstructModel(HumanGEM, samp_fullpath, genesMapping_file)
            reconstructions[individual][sample].run(media_rxns)
            # Get number of active reactions ...
            print('\n.... Get Reaction Numbers')
            for cell_type in reconstructions[individual][sample].reconstructed_models_names:
                count = 0
                for reaction in reconstructions[individual][sample].reconstructed_models[cell_type].reactions:
                    if reaction.bounds != (0, 0):
                        count += 1
                nReactions_reconstruction.loc[cell_type, sample] = count
            # Save reconstruction ...
            print('\n.... Save object and results\n')
            if not isdir(join(CRCatlasReconstruction_dir, 'NormalMatched', individual)):
                mkdir(join(CRCatlasReconstruction_dir, 'NormalMatched', individual))
            file_reconstruction_dump_object = join(CRCatlasReconstruction_dir, 'NormalMatched', individual,
                                                   ''.join(('01_', sample, '.obj')))
            with open(file_reconstruction_dump_object, 'wb') as dump_file:
                dump(reconstructions[individual][sample], dump_file)
        nReactions_reconstruction.to_csv(join(CRCatlasReconstruction_dir, 'NormalMatched', individual,
                                              '1_nReactionsReconstruction.csv'))

    '''
    Biomass GapFill
    '''

    # Optional:
    reconstructions = {}
    individuals = listdir(CRCatlasNormalMatched_dir)
    for indiv in ['SMC10']:# individuals:
        print('\n', indiv)
        reconstructions[indiv] = {}
        # Get files that starts with 02_
        indiv_samples = [file for file in listdir(join(CRCatlasReconstruction_dir, 'NormalMatched', indiv))
                         if file.startswith('01_')]
        # For each file:
        for indiv_samp in indiv_samples:
            print('- ', indiv_samp)
            samp_name = indiv_samp.replace('.obj', '').replace('01_', '')
            with open(join(CRCatlasReconstruction_dir, 'NormalMatched', indiv, indiv_samp), 'rb') as dump_file:
                temp_dump = load(dump_file)
            reconstructions[indiv][samp_name] = temp_dump

    # Gap fill:
    print('\nGap Filling models...')
    individuals = list(reconstructions.keys())
    gapfilled_models = {}
    for indiv in individuals:
        samples = list(reconstructions[indiv].keys())
        print('.. ', indiv)
        gapfilled_models[indiv] = {}

        nReactions_gafills = np.empty((len(all_cell_types), len(samples)))
        nReactions_gafills[:] = np.NaN
        nReactions_gafills = pd.DataFrame(nReactions_gafills, index=all_cell_types, columns=samples)

        for samp in samples:
            gapfilled_models[indiv][samp] = {}
            models = list(reconstructions[indiv][samp].reconstructed_models.keys())

            biomass_gapfils = np.empty((len(models), 4))
            biomass_gapfils[:] = np.NaN
            biomass_gapfils_b = pd.DataFrame(biomass_gapfils.copy(), index=models,
                                             columns=['Plasmax', 'HPLM', 'Blood_SMDB', 'Default'])
            biomass_gapfils_a = pd.DataFrame(biomass_gapfils.copy(), index=models,
                                             columns=['Plasmax', 'HPLM', 'Blood_SMDB', 'Default'])

            print('... ', samp)
            for model in models:
                print('\n..... ', model)
                gapfill_model = GapFillModel(reconstructions[indiv][samp].reconstructed_models[model],
                                             reconstructions[indiv][samp].generic_model,
                                             None, mediaOpenBounds_file)
                # Before Gapfill:
                gapfill_model.evaluate_media_biomass_capacity()
                biomass_gapfils_b.loc[model, :] = gapfill_model.media_result.iloc[:, 0]
                # Run Gapfill:
                gapfill_model.run()
                # After Gapfill:
                gapfill_model.evaluate_media_biomass_capacity()
                biomass_gapfils_a.loc[model, :] = gapfill_model.media_result.iloc[:, 0]
                count = 0
                for reaction in gapfill_model.model.reactions:
                    if reaction.bounds != (0, 0):
                        count += 1
                nReactions_gafills.loc[model, samp] = count
                gapfilled_models[indiv][samp][model] = gapfill_model.model.copy()
            biomass_gapfils_b.to_csv(join(CRCatlasReconstruction_dir, 'NormalMatched', indiv,
                                          ''.join(('3_biomassBefore_', samp, '.csv'))))
            biomass_gapfils_a.to_csv(join(CRCatlasReconstruction_dir, 'NormalMatched', indiv,
                                          ''.join(('3_biomassAfter_', samp, '.csv'))))
            file_gapfill_dump_object = join(CRCatlasReconstruction_dir, 'NormalMatched', indiv,
                                            ''.join(('02_', samp, '.obj')))
            with open(file_gapfill_dump_object, 'wb') as dump_file:
                dump(gapfilled_models[indiv][samp], dump_file)
        nReactions_gafills.to_csv(join(CRCatlasReconstruction_dir, 'NormalMatched', indiv,
                                       '2_nReactionsGapFill.csv'))

