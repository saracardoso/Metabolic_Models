if __name__ == '__main__':
    from os.path import join

    from cobra.io import read_sbml_model, write_sbml_model
    from cobra import Reaction

    from cobra.flux_analysis import find_blocked_reactions

    # Base Directories:
    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'
    models_dir = join(base_dir, '0MODELS')
    HumanGEM_dir = join(models_dir, 'HumanGEM')

    utility_data_dir = join(base_dir, 'GENERAL/utility_data')

    # --- Read HumanGEM original model ---
    print('Reading HumanGEM model...')
    HumanGEM = read_sbml_model(join(HumanGEM_dir, 'HumanGEM-1.4.1.xml.gz'))
    HumanGEM_forTcells = HumanGEM.copy()

    # --- Add the biomass reaction from the macrophage model iAB-AM0-1410:
    print('Adding biomass reaction from macrophage model iAB-AM0-1410...')
    new_reaction_mac = Reaction(id='biomass_macrophage_iABAM01410',
                                name='Biomass reaction from macrophage model iAB-AM0-1410',
                                subsystem='Artificial reactions',
                                lower_bound=0., upper_bound=1000.)
    HumanGEM_forTcells.add_reactions([new_reaction_mac])
    HumanGEM_forTcells.reactions.biomass_macrophage_iABAM01410.add_metabolites({
        'm01307c': -0.396559456, 'alpa_hs_c': -0.011499127, 'm01334c': -0.048664064, 'm01365c': -0.325724532,
        'm01369c': -0.215407845, 'm01370c': -0.282759085, 'm01371c': -25.17352552, 'm01450c': -0.020930954,
        'm01590c': -0.042373167, 'm01628c': -0.127154496,  'm00240c': -0.0036682,
        'm01639c': -0.021495345,
        'm01644c': -0.014937443, 'm01686c': -0.014937443, 'm01752c': -0.021495345, 'm01975c': -0.280436629,
        'm01974c': -0.424428935, 'm01986c': -0.366948135,
        'm03161c': -0.528027894,  # glycogen
        'm02016c': -0.043710887, 'm02040c': -25.17352552, 'm02674c': -0.004850777, 'm02675c': -0.001222285,
        'm02125c': -0.153862747, 'm02184c': -0.25953452, 'm02360c': -0.580614138, 'm02426c': -0.351852168,
        'm02471c': -0.126573882, 'm02938c': -0.004736708,  'm02646c': -0.003853116, 'm02750c': -0.003741686,
        'm02684c': -0.031527146, 'm02685c': -0.021107135, 'm02715c': -0.008918017, 'm02724c': -0.214246617,
        'm02770c': -0.346626641, 'm02808c': -0.001024655, 'm02896c': -0.476684207, 'm02908c': -0.007049706,
        'm02958c': -0.002742439,  # TAG-LD pool
        'm02993c': -0.303661194, 'm03089c': -0.069673697, 'm02494c': -0.00136164,
        'm03101c': -0.156185203, 'm03114c': -0.04602478, 'm03135c': -0.347207255,
        'm01285c': 25.17352552, 'm02039c': 25.17352552, 'm02751c': 25.17352552, 'temp001c': 1
    })
    HumanGEM_forTcells.reactions.biomass_macrophage_iABAM01410.gene_reaction_rule = ''

    # --- Remove blocked reactions
    print('Before removing blocked Reactions:')
    print('- Number of Reactions:', len(HumanGEM_forTcells.reactions))
    print('- Number of Genes:', len(HumanGEM_forTcells.genes))
    print('- Number of Metabolites:', len(HumanGEM_forTcells.metabolites))

    print('\nRemoving blocked reactions:')
    blocked_reactions = find_blocked_reactions(HumanGEM_forTcells, open_exchanges=True)
    HumanGEM_forTcells.remove_reactions(blocked_reactions, remove_orphans=True)
    print('- Number of Reactions:', len(HumanGEM_forTcells.reactions))
    print('- Number of Genes:', len(HumanGEM_forTcells.genes))
    print('- Number of Metabolites:', len(HumanGEM_forTcells.metabolites))

    # --- Write the model io a sbml file:
    print('Saving model...')
    write_sbml_model(HumanGEM_forTcells, join(HumanGEM_dir, 'HumanGEM-1.4.1_forTcells.xml.gz'))
    print('Done.')
