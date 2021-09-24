if __name__ == '__main__':
    from os.path import join
    import numpy as np
    from pandas import read_csv

    from cobra.io import read_sbml_model

    from GENERAL.code.python.ReconstructModel import ReconstructModel
    from GENERAL.code.python.EvaluateModel import GapFillModel, EvaluateModel

    # Base Directories:
    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'
    models_dir = join(base_dir, '0MODELS')
    HumanGEM_dir = join(models_dir, 'HumanGEM')

    naive_CD4_dir = join(base_dir, 'MODEL_RECONSTRUCTIONS/T_CELLS/subtypes_pipelines/naive_CD4')

    # --- Read HumanGEM_forTcells original model ---
    print('Reading HumanGEM-1.8.0...')
    HumanGEM = read_sbml_model(join(HumanGEM_dir, 'HumanGEM-1.8.0_consistent.xml.gz'))

    # Reactions in media must (?) be present:
    media_rxns = read_csv('/home/scardoso/Documents/PhD/Metabolic_Models/GENERAL/utility_data/media_openBounds.csv',
                          index_col='ID').index.tolist()

    # Get reconstructions:
    print('Reconstructing model...')
    example_rec = ReconstructModel(HumanGEM, '/home/scardoso/Documents/PhD/Metabolic_Models/0Data/scRNAseq/CRC_atlas/example.csv',
                                   '/home/scardoso/Documents/PhD/Metabolic_Models/GENERAL/utility_data/genes_mapping.json',
                                   raw_counts_cols=['Bcells', 'Tcells'])
    example_rec.run(media_rxns, '/home/scardoso/Documents/PhD/Metabolic_Models/0MODELS/example')
    # SAVE RECONSTRUCTION OBJECT (1_reconstruction)

    print('\nCharacteristics of first draft of models...')
    print('- Nº Reactions in generic model: ', len(example_rec.generic_model.reactions))
    for cell_type in example_rec.reconstructed_models_names:
        count = 0
        for reaction in example_rec.reconstructed_models[cell_type].reactions:
            if reaction.bounds != (0, 0):
                count += 1
        print('- Nº of Reactions in ', cell_type, ' draft:', count)

    # Evaluate reconstruction for biomass capacity and gapfill if necessary [gapfill_object]:
    # SAVE BEFORE AND AFTER SOLUTIONS OF BIOMASS CAPACITY IN A TABLE FOR POSTERIOR ANALYSIS TOO (2_biomass_gapfill)
    biomass_gafill = {}
    gap_filled_models = {}
    for model in example_rec.reconstructed_models.keys():
        biomass_gafill[model] = GapFillModel(example_rec.reconstructed_models[model], example_rec.generic_model,
                                             None,
                                             '/home/scardoso/Documents/PhD/Metabolic_Models/GENERAL/utility_data/media_openBounds.csv')
        # Before Gapfill:
        biomass_gafill[model].evaluate_media_biomass_capacity()
        print('\nBefore GapFill:')
        print(biomass_gafill[model].media_result)
        # Run Gapfill:
        biomass_gafill[model].run()
        # Store refinement object (in case anything needs to be evaluated again)
        # ??
        # After Gapfill:
        biomass_gafill[model].evaluate_media_biomass_capacity()
        print('\nAfter GapFill:')
        print(biomass_gafill[model].media_result)
        # Get new models, with the bounds of the added reactions as in the generic model:
        gap_filled_models[model] = biomass_gafill[model].model # SAVE THIS OBJECT (3_gapfilled_models)
        #write_sbml_model(self.reconstructed_models[col], join(out_dir, col + '.xml.gz'))

    # Check, using sampling, for expected phenotypes (organize into functions? and put them into the
    # evaluation object?) --> this will be in a separate file

