if __name__ == '__main__':
    from numpy import nan, unique
    from os.path import join
    import json
    from cobra.io import read_sbml_model

    '''
    Directories and files
    '''

    # Directories:
    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'
    HumanGEM_dir = join(base_dir, '0MODELS/HumanGEM')
    utilityData_dir = join(base_dir, 'GENERAL/utility_data')

    # Files:
    HumanGem_file = join(HumanGEM_dir, 'HumanGEM-1.8.0_consistent.xml.gz')
    reactionsSubsystems_file = join(utilityData_dir, 'reactions_subsystems.csv')

    # Subsystems IDs
    # for group in HumanGEM.groups:
    #    if group.name == 'Transport reactions': print(group.id)
    transport_group = 'group127'

    '''
    Read HumanGEM model
    '''

    print('\nReading HumanGEM-1.8.0...')
    HumanGEM = read_sbml_model(HumanGem_file)

    '''
    Create evaluations dictionary for sampling
    '''

    # 1. Calculate factor to transform media concentrations to fluxes (see create_mediumWserum.py in GENERAL/code/utils)
    print('\nCalculating factor for media transformations...')
    volume_ratio = 4000 / 176
    Tcell_weight = 60 / volume_ratio * 1e-12
    Tcell_concentration = 2.5*10e8
    Tcell_time = 48
    factor = Tcell_concentration * Tcell_weight * Tcell_time

    # 2. Create dictionary to store all different evaluation changes to test:
    print('\nCreating evaluations dictionary...')
    evals = {}
    evals['B1'] = {'change_media_bounds': {'MAR09285': 1000}, 'change_internal_bounds': None}
    evals['B2'] = {'change_media_bounds': {'MAR09045': 0}, 'change_internal_bounds': None}
    # evals['B3'] = {'change_media_bounds': {'MAR09077': 40/factor}, 'change_internal_bounds': None}  # sodium bound??
    evals['B4'] = {'change_media_bounds': {'MAR09048': 0}, 'change_internal_bounds': None}
    evals['B5'] = {'change_media_bounds': None, 'change_internal_bounds': {'MAR09517': (nan, 1000)}}
    evals['B6.1'] = {'change_media_bounds': {'MAR09135': 10 / factor}, 'change_internal_bounds': None}
    evals['B6.2'] = {'change_media_bounds': {'MAR09135': 10 / factor, 'MAR09034': '.1'}, 'change_internal_bounds': None}
    evals['B7.1'] = {'change_media_bounds': None,
                     'change_internal_bounds': {'MAR04171': (0, nan), 'MAR04193': (0, nan), 'MAR04210': (0, nan),
                                                'MAR04358': (0, nan), 'MAR04421': (0, nan), 'MAR04573': (0, nan),
                                                'MAR06627': (nan, nan), 'MAR09517': (0, nan)}}
    evals['B7.2'] = {'change_media_bounds': None,
                     'change_internal_bounds': {'MAR04171': (0, nan), 'MAR04193': (0, nan), 'MAR04210': (0, nan),
                                                'MAR04358': (0, nan), 'MAR04421': (0, nan), 'MAR04573': (0, nan),
                                                'MAR06627': (nan, nan)}}
    evals['B7.3'] = {'change_media_bounds': None,
                     'change_internal_bounds': {'MAR04171': (0, nan), 'MAR04193': (0, nan), 'MAR04210': (0, nan),
                                                'MAR04358': (0, nan), 'MAR04421': (0, nan), 'MAR04573': (0, nan),
                                                'MAR06627': (nan, nan), 'MAR09517': (nan, 1000)}}
    evals['B8.1'] = {'change_media_bounds': {'MAR09063': 1000}, 'change_internal_bounds': None}
    evals['B8.2_and_B9'] = {'change_media_bounds': {'MAR09063': 0}, 'change_internal_bounds': None}
    evals['B8.3'] = {'change_media_bounds': {'MAR09063': 1000, 'MAR09034': 0}, 'change_internal_bounds': None}
    evals['B8.4'] = {'change_media_bounds': {'MAR09063': 0, 'MAR09034': 0}, 'change_internal_bounds': None}
    evals['B8.5'] = {'change_media_bounds': {'MAR09034': 0}, 'change_internal_bounds': None}
    # --- B10:
    slc7a5_reactions = []
    slc1a5_reactions = []
    slc7a5_slc1a5_reactions = []
    for rxn in HumanGEM.groups.get_by_id(transport_group).members:
        for gene in rxn.genes:
            if gene.id == 'ENSG00000103257':
                slc7a5_reactions = slc7a5_reactions + [rxn.id]
            elif gene.id == 'ENSG00000103257':
                slc1a5_reactions = slc1a5_reactions + [rxn.id]
    slc7a5_slc1a5_reactions = list(unique(slc7a5_reactions + slc1a5_reactions))
    evals['B10.1'] = {'change_media_bounds': None, 'change_internal_bounds': {}}
    for rxn in slc7a5_reactions:
        evals['B10.1']['change_internal_bounds'][rxn] = (0, 0)
    evals['B10.2'] = {'change_media_bounds': None, 'change_internal_bounds': {}}
    for rxn in slc1a5_reactions:
        evals['B10.2']['change_internal_bounds'][rxn] = (0, 0)
    evals['B10.3'] = {'change_media_bounds': None, 'change_internal_bounds': {}}
    for rxn in slc7a5_slc1a5_reactions:
        evals['B10.3']['change_internal_bounds'][rxn] = (0, 0)
    # ---
    evals['B11'] = {'change_media_bounds': {'MAR09069': 0, 'MAR09067': 0}, 'change_internal_bounds': None}
    evals['B12.1'] = {'change_media_bounds': {'MAR09066': 0}, 'change_internal_bounds': None}
    evals['B12.2'] = {'change_media_bounds': {'MAR09065': 0}, 'change_internal_bounds': None}
    evals['B12.3'] = {'change_media_bounds': {'MAR09066': 0, 'MAR09065': 0}, 'change_internal_bounds': None}
    evals['B13'] = {'change_media_bounds': None, 'change_internal_bounds': {'MAR06916': (0, 0)}}
    evals['B13.o'] = {'change_media_bounds': {'MAR09048': 0}, 'change_internal_bounds': {'MAR06916': (0, 0)}}
    evals['B14'] = {'change_media_bounds': None, 'change_internal_bounds': {'MAR09517': (0, 0)}}
    evals['B14.o'] = {'change_media_bounds': {'MAR09048': 0}, 'change_internal_bounds': {'MAR09517': (0, 0)}}
    evals['B15.1'] = {'change_media_bounds': None, 'change_internal_bounds': {'MAR00030': (0, 0), 'MAR04156': (0, 0),
                                                                              'MAR07673': (0, 0)}}
    evals['B15.2'] = {'change_media_bounds': None, 'change_internal_bounds': {'MAR04156': (0, 0), 'MAR04295': (0, 0),
                                                                              'MAR07673': (0, 0)}}
    evals['B15.3'] = {'change_media_bounds': None, 'change_internal_bounds': {'MAR00030': (0, 0), 'MAR04156': (0, 0),
                                                                              'MAR04295': (0, 0), 'MAR07673': (0, 0)}}
    evals['B16'] = {'change_media_bounds': {'MAR09048': 0}, 'change_internal_bounds': None}
    evals['B17'] = {'change_media_bounds': {'MAR09262': 0, 'MAR09275': 0, 'MAR09220': 0, 'MAR09436': 0, 'MAR09343': 0,
                                            'MAR01923': 0, 'MAR09847': 0}, 'change_internal_bounds': None}
    evals['B21'] = {'change_media_bounds': None, 'change_internal_bounds': {'MAR03839': (0, 0)}}

    # 3. Save dictionary in json file:
    print('\nSaving evaluations into a json file...')
    with open(join(utilityData_dir, 'sampling_evals_Tcells.json'), 'w') as f:
        json.dump(evals, f)

    print('\nDone!')
