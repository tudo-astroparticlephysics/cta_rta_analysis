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
from tqdm import tqdm
from collections import Counter
from itertools import chain
import os


# do some horrible things to silencece astropy warnings in ctapipe
warnings.filterwarnings('ignore', category=AstropyDeprecationWarning, append=True)
warnings.filterwarnings('ignore', category=FutureWarning, append=True)

SubMomentParameters = namedtuple('SubMomentParameters', 'size,cen_x,cen_y,length,width,psi')


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
@click.option('-n', '--n_jobs', help='Number of threads to use', default=-1)
@click.option('-t', '--tel_type', help='Telescope Types to use', type=click.Choice(['all', 'MST', 'SST', 'LST']), default='all')
@click.option('-y', '--yes', help='Override all prompts. Overwrites exisitng files', default=False, is_flag=True)
@click.option('-l', '--event_limit', help='Number of events to reconstruct', default=-1,)
def main(input_file_path, output_file_path, instrument_description, n_jobs, tel_type, yes, event_limit):


    if os.path.exists(output_file_path):
        if not yes:
            click.confirm(f'File {output_file_path} exists. Overwrite?', default=False, abort=True)
        os.remove(output_file_path)

    instrument = pickle.load(open(instrument_description, 'rb'))

    runs = fact.io.read_data(input_file_path, key='runs')

    telescope_events = fact.io.read_data(input_file_path, key='telescope_events',)
    telescope_events.set_index(['run_id', 'array_event_id', 'telescope_id'], inplace=True, drop=False, verify_integrity=True)

    if tel_type != 'all':
        telescope_events = telescope_events.query(f'telescope_type_name == "{tel_type}"')

    array_events = fact.io.read_data(input_file_path, key='array_events',)
    array_events.set_index(['run_id', 'array_event_id'], verify_integrity=True, inplace=True, drop=False)

    events = pd.merge(left=array_events, right=telescope_events, left_on=['run_id', 'array_event_id'], right_on=['run_id', 'array_event_id']).dropna()
    if event_limit > 1:
        events = events[0:event_limit]
        print(f'Reconstruncting {len(events) events.}')

    if n_jobs == -1:
        n_jobs = multiprocessing.cpu_count() // 2

    if n_jobs > 1:
        results = Parallel(n_jobs=n_jobs, verbose=5)(delayed(reconstruct_direction)(array_event_id, group, instrument=instrument) for array_event_id, group in events.groupby('array_event_id'))
    else:
        results = [reconstruct_direction(array_event_id, group, instrument=instrument) for array_event_id, group in tqdm(events.groupby(['run_id', 'array_event_id']))]


    telescope_features = pd.DataFrame(list(chain.from_iterable([r[1] for r in results if r is not None])))
    telescope_features.set_index(['run_id', 'array_event_id', 'telescope_id'], inplace=True)

    array_features = pd.DataFrame([r[0] for r in results if r is not None])
    array_features.set_index(['run_id', 'array_event_id'], inplace=True)

    columns_to_delete = set(array_features.columns) & set(array_events.columns)
    array_events.drop(columns=columns_to_delete, inplace=True)
    array_events = pd.merge(left=array_features, right=array_events, left_index=True, right_index=True)

    columns_to_delete = set(telescope_features.columns) & set(telescope_events.columns)
    telescope_events.drop(columns=columns_to_delete, inplace=True)
    telescope_events = pd.merge(left=telescope_features, right=telescope_events, left_index=True, right_index=True)

    if 'gamma_prediction' in telescope_events.columns:
        array_events['gamma_prediction_mean'] = telescope_events.groupby(['run_id', 'array_event_id'])['gamma_prediction'].mean()
        array_events['gamma_prediction_std'] = telescope_events.groupby(['run_id', 'array_event_id'])['gamma_prediction'].std()
    if 'gamma_energy_prediction' in telescope_events.columns:
        array_events['gamma_energy_prediction_mean'] = telescope_events.groupby(['run_id', 'array_event_id'])['gamma_energy_prediction'].mean()
        array_events['gamma_energy_prediction_std'] = telescope_events.groupby(['run_id', 'array_event_id'])['gamma_energy_prediction'].std()


    fact.io.write_data(runs, output_file_path, key='runs')
    array_events.reset_index(drop=True, inplace=True)
    fact.io.write_data(array_events, output_file_path, key='array_events', mode='a')
    telescope_events.reset_index(drop=True, inplace=True)
    fact.io.write_data(telescope_events, output_file_path, key='telescope_events', mode='a')


def reconstruct_direction(array_event_id, group, instrument):

    reco = HillasReconstructor()


    params = {}
    pointing_azimuth = {}
    pointing_altitude = {}
    for index, row in group.iterrows():
        tel_id = row.telescope_id
        # the data in each event has to be put inside these namedtuples to call reco.predict
        moments = SubMomentParameters(size=row.intensity, cen_x=row.x * u.m, cen_y=row.y * u.m, length=row.length * u.m, width=row.width * u.m, psi=row.psi * u.rad)
        params[tel_id] = moments
        pointing_azimuth[tel_id] = row.pointing_azimuth * u.rad
        pointing_altitude[tel_id] = ((np.pi/2) - row.pointing_altitude) * u.rad
# ((np.pi/2) - event.mc.tel[telescope_id].altitude_raw )* u.rad
    try:
        reconstruction = reco.predict(params, instrument, pointing_azimuth, pointing_altitude)
    except NameError:
        return None

    c = Counter(group.telescope_type_name)
    total_signal = group.intensity.sum()

    event_features = {'alt_prediction': reconstruction.alt.si.value,  # TODO srsly now? FFS
            'az_prediction': reconstruction.az.si.value,
            'core_x_prediction': reconstruction.core_x.si.value,
            'core_y_prediction': reconstruction.core_y.si.value,
            'array_event_id': array_event_id,
            'run_id': group.run_id.iloc[0],
            'h_max_prediction': reconstruction.h_max.si.value,
            'total_intensity': total_signal,
            'num_triggered_lst': c['LST'],
            'num_triggered_mst': c['MST'],
            'num_triggered_sst': c['SST'],
            }

    tel_wise = []
    for index, row in group.iterrows():
        tel_id = row.telescope_id
        f = {}
        pos = instrument.subarray.positions[tel_id]
        x, y = pos[0], pos[1]
        core_x = reconstruction.core_x
        core_y = reconstruction.core_y
        d = np.sqrt((core_x - x)**2 + (core_y - y)**2)
        f['distance_to_core'] = d.si.value

        t = instrument.subarray.tel[tel_id]
        f['area'] = t.optics.mirror_area.si.value
        f['focal_length'] = t.optics.equivalent_focal_length.si.value

        f['normed_intensity'] = row.intensity / t.optics.mirror_area.si.value
        f['scaled_width'] = row.width / t.optics.equivalent_focal_length.si.value
        f['array_event_id'] = array_event_id
        f['telescope_id'] = tel_id
        f['run_id'] = row.run_id
        tel_wise.append(f)


    return event_features, tel_wise


if __name__ == '__main__':
    main()
