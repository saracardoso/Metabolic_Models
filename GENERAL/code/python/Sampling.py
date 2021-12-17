from pandas import read_csv
from datetime import datetime

from cobra.sampling import OptGPSampler, ACHRSampler


def sampling(model, media_file, change_media_bounds: dict = None, change_internal_bounds: dict = None,
             n_samples: int = 1000, processes: int = 4):
    # if number in a change_media_bounds' value is a string, it will be multiplied by the value in the medium
    media = read_csv(media_file, index_col='ID')
    sampling_results = {}
    with model as model_test:
        for medium in media.columns:
            print('\n', datetime.now().strftime('%H:%M:%S'), ' | ', 'MEDIUM: ', medium)
            model_test.medium = media[medium].to_dict()
            if change_media_bounds is not None:
                print(datetime.now().strftime('%H:%M:%S'), ' | ', '- Changing media bounds..')
                new_media = model_test.medium.copy()
                for cpd_in, val in change_media_bounds[medium].items():
                    if isinstance(val, str): new_media[cpd_in] = float(val) * new_media[cpd_in]
                    else: new_media[cpd_in] = val
                model_test.medium = new_media.copy()
            if change_internal_bounds is not None:
                print(datetime.now().strftime('%H:%M:%S'), ' | ', '- Changing internal bounds..')
                for rxn, rxn_bounds in change_internal_bounds.values():
                    model_test.reactions.get_by_id(rxn).bounds = rxn_bounds
            print(datetime.now().strftime('%H:%M:%S'), ' | ', 'Initializing sampler...')
            sampler = OptGPSampler(model, processes=processes, thinning=100)
            # sampler = ACHRSampler(model, thinning=100)
            print(datetime.now().strftime('%H:%M:%S'), ' | ', 'Sampling...')
            try:
                fluxes_samples = sampler.sample(n_samples)
            except RuntimeError:
                sampling_results[medium] = None
                print('Did not work (Can not escape sampling region, model seems numerically unstable).')
            else:
                print(datetime.now().strftime('%H:%M:%S'), ' | ', 'Getting validation status...')
                val_status = sampler.validate(fluxes_samples)
                print(datetime.now().strftime('%H:%M:%S'), ' | ', 'Done!')
                print(datetime.now().strftime('%H:%M:%S'), ' | ', 'Filtering invalid samples...')
                sampling_results[medium] = fluxes_samples[val_status == 'v']
            print(datetime.now().strftime('%H:%M:%S'), ' | ', 'Done!')
    return sampling_results

