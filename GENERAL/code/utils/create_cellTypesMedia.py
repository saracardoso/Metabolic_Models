if __name__ == '__main__':
    from os.path import join

    from GENERAL.code.python.others import construct_media_with_serum, construct_media_fluxes_file

    '''
    Directories and files
    '''

    # Directories:
    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'
    utilityData_dir = join(base_dir, 'GENERAL/utility_data')

    # Files:
    media_file = join(utilityData_dir, 'media.csv')

    '''
    Create media with fluxes specific for cell-types

    Convertion of media concentrations to fluxes is performed using the following equation:
        flux = ([metabolite]]) / ([cell] * cell_weight * time)

        [metabolite] --> concentration of the metabolite in the medium, in mM
        [cell]       --> concentration of cells, in nºcells/L
        cell_weight  --> dry weight of a cell, in gDW
        time         --> time that the cell takes to consume the metabolite
    '''

    '''
    a. Tcells
    '''

    # Knowing that osteosarcoma (U2OS) cells have a cell dry weight of ~60 pg and a cell volume of 4000 um3
    # And that the volume of a T lymphocyte is ~176 um3
    # The T lymphocyte dry weight can be calculated as follows:
    volume_ratio = 4000 / 176  # which makes 22.73 ratio
    Tcell_weight = 60 / volume_ratio * 1e-12  # which makes ~2.640*10e-12 gDW

    # Considering the optimal Tcell concentration in a medium to be 2.5*10e5 cells/mL
    Tcell_concentration = 2.5 * 10e8  # cells/L

    # Considering that culture media is normally changed every two days:
    Tcell_time = 48  # h

    # The following metabolites will be set to unconstrained, i.e, always available:
    unconstrain = ['MAR09047',  # H2O
                   'MAR09048',  # O2
                   'MAR09079'  # H+
                   ]

    # Media concentration values will be converted to fluxes using these values and saved to a file:
    # Plus, columns will be added to account for all the media metabolites to be unconstrained.
    construct_media_fluxes_file(join(utilityData_dir, 'media_8percSerum.csv'),
                                join(utilityData_dir, 'media_Tcells_8percSerum.csv'),
                                Tcell_concentration, Tcell_weight, Tcell_time,
                                unconstrained_metabolites=unconstrain, add_open_bounds_column=False,
                                media_to_use=['Plasmax_serum', 'Blood_SMDB'])

    '''
    b. CRC cells
    '''
