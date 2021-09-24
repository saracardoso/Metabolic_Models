if __name__ == '__main__':
    from os.path import join
    from os import getcwd

    from GENERAL.code.python.others import construct_media_fluxes_file

    # Base Directories:
    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'

    # Convertion of media concentrations to fluxes is performed using the following equation:
    # flux = ([metabolite]]) / ([cell] * cell_weight * time)
    #
    # [metabolite] --> concentration of the metabolite in the medium, in mM
    # [cell]       --> concentration of cells, in nºcells/L
    # cell_weight  --> dry weight of a cell, in gDW
    # time         --> time that the cell takes to consume the metabolite
    #
    # Values calculated elsewhere

    # The average dry weight of a CRC cell is:
    CRC_weight = 3.96632116942834*10e-11  # gDW

    # Considering the optimal CRC concentration in a medium to be 1*10e7 cells/L
    CRC_concentration = 1.35*10e7  # cells/L

    # Considering that culture media is normally changed every two days:
    CRC_time = 48  # h

    # The following metabolites will be set to unconstrained, i.e, always available:
    unconstrain = ['HMR_9047',  # H2O
                   'HMR_9048',  # O2
                   'HMR_9079'  # H+
                   ]

    # Media concentration values will be converted to fluxes using these values and saved to a file:
    # Plus, columns will be added to account for all the media metabolites to be unconstrained.
    CRC_media_fluxes = join(base_dir, 'MODEL_RECONSTRUCTIONS/CRC/general/utility_data/CRC_media_fluxes.csv')
    construct_media_fluxes_file(base_dir, CRC_media_fluxes, CRC_concentration, CRC_weight, CRC_time,
                                unconstrained_metabolites=unconstrain, add_open_bounds_column=True,
                                media_to_use=['Plasmax', 'Blood_SMDB'])
