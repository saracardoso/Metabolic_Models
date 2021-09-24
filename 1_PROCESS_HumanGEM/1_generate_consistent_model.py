if __name__ == '__main__':
    from os.path import join

    from cobra.io import read_sbml_model, write_sbml_model

    from cobra.flux_analysis import find_blocked_reactions

    # Base Directories:
    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'
    models_dir = join(base_dir, '0MODELS')
    HumanGEM_dir = join(models_dir, 'HumanGEM')

    utility_data_dir = join(base_dir, 'GENERAL/utility_data')
    test_humanGEM_utility_data_dir = join(base_dir, '1_PROCESS_HumanGEM/utility_data')

    # --- Read HumanGEM original model ---
    print('Reading HumanGEM-1.8.0:')
    HumanGEM = read_sbml_model(join(HumanGEM_dir, 'HumanGEM-1.8.0.xml.gz'))
    print('- Number of Reactions:', len(HumanGEM.reactions))
    print('- Number of Genes:', len(HumanGEM.genes))
    print('- Number of Metabolites:', len(HumanGEM.metabolites))

    # --- Remove blocked reactions ---
    print('\nRemoving blocked reactions:')
    blocked_reactions = find_blocked_reactions(HumanGEM, open_exchanges=True)
    HumanGEM.remove_reactions(blocked_reactions, remove_orphans=True)
    print('- Number of Reactions:', len(HumanGEM.reactions))
    print('- Number of Genes:', len(HumanGEM.genes))
    print('- Number of Metabolites:', len(HumanGEM.metabolites))

    # --- Saving consistent model:
    print('\nSaving consistent HumanGEM model...')
    write_sbml_model(HumanGEM, join(HumanGEM_dir, 'HumanGEM-1.8.0_consistent.xml.gz'))

    print('\nDone.')
