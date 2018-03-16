import click
from ctapipe.reco import HillasReconstructor
import pickle
import fact.io
from collections import namedtuple
from tqdm import tqdm
import astropy.units as u
import pandas as pd
import warnings
from astropy.utils.exceptions import AstropyDeprecationWarning


warnings.filterwarnings('ignore', category=AstropyDeprecationWarning, append=True)
warnings.filterwarnings('ignore', category=FutureWarning, append=True)

SubMomentParameters = namedtuple('SubMomentParameters', 'size,cen_x,cen_y,length,width,psi')


@click.command()
@click.argument(
    'predicted_events', type=click.Path(
        exists=True,
        dir_okay=False,
    ))
@click.argument(
    'instrument_description', type=click.Path(
        exists=True,
        dir_okay=False,
    ))
@click.argument(
    'output_file', type=click.Path(
        exists=False,
        dir_okay=False,
    ))
def main(predicted_events, instrument_description, output_file):
    reco = HillasReconstructor()
    instrument = pickle.load(open(instrument_description, 'rb'))
    events = fact.io.read_data(predicted_events, key='events')

    results = []
    for array_event_id, group in tqdm(events.groupby('array_event_id')):
        params = {}
        pointing_azimuth = {}
        pointing_altitude = {}
        for index, row in group.iterrows():
            tel_id = row.telescope_id
            moments = moments_from_row(row)
            params[tel_id] = moments
            pointing_azimuth[tel_id] = row.pointing_azimuth * u.rad
            pointing_altitude[tel_id] = row.pointing_altitude * u.rad

        reconstruction = reco.predict(params, instrument, pointing_azimuth, pointing_altitude)
        results.append(create_result(array_event_id, group, reconstruction))

    df = pd.DataFrame(results)
    fact.io.write_data(df, output_file, key='events', mode='w')


def create_result(array_event_id, group, reconstruction):
    d = {
        'alt_prediction': reconstruction.alt.si.value,
        'az_prediction': reconstruction.az.si.value,
        'gamma_energy_prediction': group.gamma_energy_prediction.mean(),
        'gamma_prediction': group.gamma_prediction.mean(),
        'core_x_prediction': reconstruction.core_x.si.value,
        'core_y_prediction': reconstruction.core_y.si.value,
        'h_max_prediction': reconstruction.h_max.si.value,
        'array_event_id': array_event_id,
        'mc_alt': group.mc_alt.mean(),
        'mc_az': group.mc_az.mean(),
        'mc_core_x': group.mc_core_y.mean(),
        'mc_core_y': group.mc_core_x.mean(),
        'mc_energy': group.mc_energy.mean(),
        'mc_height_first_interaction': group.mc_height_first_interaction.mean(),
        'mc_corsika_primary_id': int(group.mc_corsika_primary_id.mean()),
        'pointing_azimuth': group.pointing_azimuth.mean(),
        'pointing_altitude': group.pointing_altitude.mean(),
    }
    return d


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
