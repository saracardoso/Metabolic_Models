if __name__ == '__main__':
    from os.path import join
    from os import listdir
    from pandas import read_csv

    from cobra.io import read_sbml_model

    from GENERAL.code.python.ReconstructModel import ReconstructModel
    from GENERAL.code.python.EvaluateModel import GapFillModel

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

    print('Reading HumanGEM-1.8.0...')
    HumanGEM = read_sbml_model(HumanGem_file)

    '''
    Get reactions present in media
    '''

    print('Getting media reactions...')
    media_rxns = read_csv(mediaOpenBounds_file, index_col='ID').index.tolist()

    '''
    Reconstruct model
    '''

    print('Reconstructing models...')
    individuals = listdir(CRCatlasNormalMatched_dir)
    reconstructions = {}
    for individual in individuals:
        reconstructions[individual] = {}
        samples = listdir(join(CRCatlasNormalMatched_dir, individual))
        for samp in samples:
            samp_fullpath = join(CRCatlasNormalMatched_dir, individual, samp)
            samp_name = samp.replace('.csv', '')
            reconstructions[individual][samp] = ReconstructModel(HumanGEM, samp_fullpath, genesMapping_file)
            reconstructions[individual][samp].run(media_rxns)
            # Save reconstruction ...

    # [CONTINUE...]