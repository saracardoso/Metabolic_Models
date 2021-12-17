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
    Create media file considering serum metabolites at 8%
    '''

    print('Creating media with 8% serum...')
    construct_media_with_serum(media_file, join(utilityData_dir, 'media_8percSerum.csv'), serum_percentage=8,
                               serum_medium='Blood_SMDB', media_to_change=['Plasmax', 'HPLM'])

    print('Done!')
