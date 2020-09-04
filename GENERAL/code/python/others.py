from pandas import read_csv
from os.path import join
import warnings

import numpy as np


def construct_media_fluxes_file(metabolic_models_dir, path_outFile, cell_concentration, cell_weight, time,
                                unconstrained_metabolites=None, add_open_bounds_column=True,
                                media_to_use=['Plasmax', 'HPLM', 'RPMI-1640', 'Blood_SMDB']):
    media_concentrations = read_csv(join(metabolic_models_dir, 'GENERAL/utility_data/media.csv'), index_col='ID')

    # Check if all media fiven in media_to_use exist in file:
    not_present = []
    for medium_name in media_to_use:
        if medium_name not in media_concentrations.columns:
            not_present = not_present + [medium_name]
    if len(not_present) > 0:
        warnings.warn('The following media names are invalid: ', not_present)
        return
    # print('All media names given are valid.')

    # Convert metabolites present in the media to fluxes:
    media_fluxes = media_concentrations.loc[:, media_to_use]
    factor = cell_concentration * cell_weight * time
    media_fluxes = media_fluxes / factor

    # Add unconstrained metabolites:
    if unconstrained_metabolites is not None:
        for metab in unconstrained_metabolites:
            nrow, ncol = media_fluxes.shape
            media_fluxes.loc[metab] = [1000] * ncol

    # Add columns for media without the concentrations constraint:
    if add_open_bounds_column:
        for medium in media_to_use:
            to_add = np.array([0] * media_fluxes.shape[0])
            where_open = np.where(media_fluxes.loc[:, medium] != 0)[0]
            to_add[where_open] = 1000
            new_name = medium + '_open'
            media_fluxes.loc[:, new_name] = to_add.tolist()

    # Save final fluxes into the file:
    media_fluxes.to_csv(path_or_buf=path_outFile)
