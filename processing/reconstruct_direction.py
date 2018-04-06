import click
from ctapipe.reco import HillasReconstructor
import pickle
import fact.io
from collections import namedtuple
import astropy.units as u
import pandas as pd
import warnings
from astropy.utils.exceptions import AstropyDeprecationWarning
from joblib import Parallel, delayed
import numpy as np
import multiprocessing


# do some horrible things to silencece astropy warnings in ctapipe
warnings.filterwarnings('ignore', category=AstropyDeprecationWarning, append=True)
warnings.filterwarnings('ignore', category=FutureWarning, append=True)

SubMomentParameters = namedtuple('SubMomentParameters', 'size,cen_x,cen_y,length,width,psi')


def dummy_function_h_max(self, hillas_dict, subarray, tel_phi):
    return -1


@click.command()
@click.argument(
    'input_file_path', type=click.Path(
        exists=True,
        dir_okay=False,
    ))
@click.argument(
    'output_file_path', type=click.Path(
        exists=False,
        dir_okay=False,
    ))
@click.argument(
    'instrument_description', type=click.Path(
        exists=True,
        dir_okay=False,
    ))
@click.option('-y', '--yes', help='Do not prompt for overwrites', is_flag=True)
def main(input_file_path, output_file_path, instrument_description, yes):

    instrument = pickle.load(open(instrument_description, 'rb'))

    runs = fact.io.read_data(input_file_path, key='runs')
    runs.set_index('run_id', drop=True, verify_integrity=True, inplace=True)

    telescope_events = fact.io.read_data(input_file_path, key='telescope_events')
    telescope_events.set_index('telescope_event_id', drop=True, verify_integrity=True, inplace=True)

    array_events = fact.io.read_data(input_file_path, key='array_events')
    array_events.set_index('array_event_id', drop=True, verify_integrity=True, inplace=True)


    events = pd.merge(left=array_events, right=telescope_events, left_index=True, right_on='array_event_id')

    n_jobs = multiprocessing.cpu_count()
    results = Parallel(n_jobs=n_jobs, verbose=5) (delayed(reconstruct_direction)(array_event_id, group, instrument=instrument) for array_event_id, group in events.groupby('array_event_id'))
    assert len(results) == len(array_events)
    df = pd.DataFrame(results)
    df.set_index('array_event_id', inplace=True)

    array_events['alt_prediction'] = df.alt_prediction
    array_events['az_prediction'] = df.az_prediction
    array_events['core_x_prediction'] = df.core_x_prediction
    array_events['core_y_prediction'] = df.core_y_prediction
    if 'gamma_prediction' in telescope_events.columns:
        array_events['gamma_prediction_mean'] = telescope_events.groupby('array_event_id')['gamma_prediction'].mean()
        array_events['gamma_prediction_std'] = telescope_events.groupby('array_event_id')['gamma_prediction'].std()
    if 'gamma_energy_prediction' in telescope_events.columns:
        array_events['gamma_energy_prediction_mean'] = telescope_events.groupby('array_event_id')['gamma_energy_prediction'].mean()
        array_events['gamma_energy_prediction_std'] = telescope_events.groupby('array_event_id')['gamma_energy_prediction'].std()


    fact.io.write_data(runs, output_file_path, key='runs')
    fact.io.write_data(array_events, output_file_path, key='array_events', mode='a')
    fact.io.write_data(telescope_events, output_file_path, key='telescope_events', mode='a')


def reconstruct_direction(array_event_id, group, instrument):

    reco = HillasReconstructor()
    # monkey patch this huansohn. this is super slow otherwise. who needs max h anyways
    reco.fit_h_max = dummy_function_h_max

    params = {}
    pointing_azimuth = {}
    pointing_altitude = {}
    for index, row in group.iterrows():
        tel_id = row.telescope_id
        # the data in each event has to be put inside these namedtuples to call reco.predict
        moments = SubMomentParameters(size=row.intensity, cen_x=row.x * u.m, cen_y=row.y * u.m, length=row.length * u.m, width=row.width * u.m, psi=row.psi * u.rad)
        params[tel_id] = moments
        pointing_azimuth[tel_id] = row.pointing_azimuth * u.rad
        pointing_altitude[tel_id] = row.pointing_altitude * u.rad


    reconstruction = reco.predict(params, instrument, pointing_azimuth, pointing_altitude)
    if reconstruction.alt.si.value == np.nan:
        print('Not reconstructed')
        print(params)

    return {'alt_prediction': ((np.pi / 2) - reconstruction.alt.si.value),  # TODO srsly now? FFS
            'az_prediction': reconstruction.az.si.value,
            'core_x_prediction': reconstruction.core_x.si.value,
            'core_y_prediction': reconstruction.core_y.si.value,
            'array_event_id': array_event_id,
            # 'h_max_prediction': reconstruction.h_max.si.value
            }


def strip_unit(v):
    try:
        return v.si.value
    except AttributeError:
        return v


def moments_from_row(row):
    return SubMomentParameters(
        size=row.intensity,
        cen_x=row.x * u.m,
        cen_y=row.y * u.m,
        length=row.length * u.m,
        width=row.width * u.m,
        psi=row.psi * u.rad
    )


if __name__ == '__main__':
    main()
