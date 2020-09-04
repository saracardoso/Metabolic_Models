if __name__ == '__main__':
    from os.path import join

    from cobra.io import read_sbml_model

    from GENERAL.code.python.EvaluateModel import EvaluateModel

    # Base Directories:
    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'
    models_dir = join(base_dir, '0MODELS')
    HumanGEM_dir = join(models_dir, 'HumanGEM')

    CRC_utility_data_dir = join(base_dir, 'MODEL_RECONSTRUCTIONS/CRC/general/utility_data')

    # Files:
    media_file = join(CRC_utility_data_dir, 'CRC_media_fluxes.csv')

    # --- Read HumanGEM_forTcells original model ---
    print('Reading HumanGEM-1.4.1...')
    HumanGEM = read_sbml_model(join(HumanGEM_dir, 'HumanGEM-1.4.1.xml.gz'))

    # --- Evaluate Model ---
    evaluate_HumanGEM_forCRC = EvaluateModel(model=HumanGEM, media_file=media_file)
    # - Media biomass capacity:
    print('Evaluating model capacity to produce each biomass reaction from the different media...')
    evaluate_HumanGEM_forCRC.evaluate_media_biomass_capacity()
    evaluate_HumanGEM_forCRC.save_media_biomass_capacity_csv(
        join(base_dir, 'MODEL_RECONSTRUCTIONS/CRC/media_biomass_capacity_results.csv'))
    print('Done.')
