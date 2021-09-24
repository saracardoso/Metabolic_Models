if __name__ == '__main__':
    from os.path import join

    from cobra.io import read_sbml_model

    from GENERAL.code.python.EvaluateModel import EvaluateModel

    # Base Directories:
    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'
    models_dir = join(base_dir, '0MODELS')
    HumanGEM_dir = join(models_dir, 'HumanGEM')

    utility_data_dir = join(base_dir, 'GENERAL/utility_data')

    # Files:
    metabolic_tasks_hsa_all = join(utility_data_dir, 'metabolic_tasks_hsa_all_consistent.json')
    media_file = join(utility_data_dir, 'media_openBounds.csv')

    # --- Read HumanGEM original model ---
    print('Reading HumanGEM-1.8.0...')
    HumanGEM = read_sbml_model(join(HumanGEM_dir, 'HumanGEM-1.8.0_consistent.xml.gz'))

    # --- Evaluate Model ---
    print('\nInitializing model evaluation:')
    evaluate_HumanGEM = EvaluateModel(model=HumanGEM, tasks_file=metabolic_tasks_hsa_all,
                                      media_file=media_file)
    # - Task evaluation:
    print('Evaluating model capacity to perform all the metabolic tasks...')
    evaluate_HumanGEM.evaluate_tasks()
    evaluate_HumanGEM.save_tasks_result_csv(
        join(base_dir, '1_PROCESS_HumanGEM/results/1_tasks_consistent_model.csv'))
    # - Media biomass capacity:
    print('Evaluating model capacity to produce each biomass reaction from the different media...')
    evaluate_HumanGEM.evaluate_media_biomass_capacity()
    evaluate_HumanGEM.save_media_biomass_capacity_csv(
        join(base_dir, '1_PROCESS_HumanGEM/results/1_media_consistent_model.csv'))

    print('\nDone.')