if __name__ == '__main__':
    from os.path import join

    from cobra.io import read_sbml_model

    from GENERAL.code.python.EvaluateModel import EvaluateModel

    # Base Directories:
    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'
    models_dir = join(base_dir, '0MODELS')
    HumanGEM_dir = join(models_dir, 'HumanGEM')

    Tcells_utility_data_dir = join(base_dir, 'MODEL_RECONSTRUCTIONS/T_CELLS/general/utility_data')

    # Files:
    media_file = join(Tcells_utility_data_dir, 'Tcell_media_fluxes.csv')

    # --- Read HumanGEM_forTcells original model ---
    print('Reading HumanGEM-1.4.1_forTcells...')
    HumanGEM_forTcells = read_sbml_model(join(HumanGEM_dir, 'HumanGEM-1.4.1_forTcells.xml.gz'))

    # --- Evaluate Model ---
    evaluate_HumanGEM_forTcells = EvaluateModel(model=HumanGEM_forTcells, media_file=media_file)
    # - Media biomass capacity:
    print('Evaluating model capacity to produce each biomass reaction from the different media...')
    evaluate_HumanGEM_forTcells.evaluate_media_biomass_capacity()
    evaluate_HumanGEM_forTcells.save_media_biomass_capacity_csv(
        join(base_dir, 'MODEL_RECONSTRUCTIONS/T_CELLS/media_biomass_capacity_results.csv'))
    print('Done.')
