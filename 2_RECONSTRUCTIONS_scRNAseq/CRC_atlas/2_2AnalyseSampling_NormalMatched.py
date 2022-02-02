if __name__ == '__main__':
    from os.path import join
    from os import listdir
    from dill import load
    import json
    from pandas import concat, Series
    from gc import collect
    from re import search

    from cobra.io import read_sbml_model

    '''
    Directories and files
    '''

    # Directories:
    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'
    CRCReconstruction_dir = join(base_dir, '2_RECONSTRUCTIONS_scRNAseq/CRC_atlas')
    CRCReconstructionNormalMatched_dir = join(CRCReconstruction_dir, 'NormalMatched')
    HumanGEM_dir = join(base_dir, '0MODELS/HumanGEM')

    # Files:
    CRCatlas_sampling_file = join(CRCReconstruction_dir, 'CRCatlas_sampling.json')
    HumanGem_file = join(HumanGEM_dir, 'HumanGEM-1.8.0_consistent.xml.gz')

    '''
    Read HumanGEM model
    '''

    print('\nReading HumanGEM-1.8.0...')
    HumanGEM = read_sbml_model(HumanGem_file)

    '''
    Get number of reactions in the models
    '''

    CRCatlas_sampling = json.load(open(CRCatlas_sampling_file))
    individuals = listdir(CRCReconstructionNormalMatched_dir)
    reactions_df = None
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
                    print('--', cell_type)
                    if cell_type in CRCatlas_sampling['NormalMatched'][indiv]['control'][samp_name]:
                        rxn_vec = [0] * len(model.reactions)
                        for idx, rxn in enumerate(model.reactions):
                            if rxn.bounds != (0, 0):
                                rxn_vec[idx] = 1
                        models_names.append('_'.join((indiv, samp_name, cell_type)))
                        if reactions_df is None:
                            reactions_df = concat([Series(rxn_vec)], axis=1)
                        else:
                            reactions_df = concat([reactions_df, Series(rxn_vec)], axis=1)
                del temp_dump
                collect()
    reactions_df.columns = models_names
    reactions_df.to_csv(join(CRCReconstructionNormalMatched_dir, '1_control_analysis/reaction_presence.csv'))

