from ctapipe.io.eventsourcefactory import EventSourceFactory
from ctapipe.calib import CameraCalibrator
from ctapipe.image.hillas import hillas_parameters
from ctapipe.image.cleaning import tailcuts_clean
import pandas as pd
import fact.io
import click
import os
import pyhessio
import numpy as np

names_to_id = {'LSTCam': 1, 'NectarCam': 2, 'FlashCam': 3, 'DigiCam': 4, 'CHEC': 5}
types_to_id = {'LST': 1, 'MST': 2, 'SST': 3}
allowed_cameras = ['LSTCam', 'NectarCam', 'DigiCam']


cleaning_level = {
                    'ASTRICam': (5, 7),  # (5, 10)?
                    'FlashCam': (12, 15),
                    'LSTCam': (5, 10),  # ?? (3, 6) for Abelardo...
                    # ASWG Zeuthen talk by Abelardo Moralejo:
                    'NectarCam': (4, 8),
                    # "FlashCam": (4, 8),  # there is some scaling missing?
                    'DigiCam': (3, 6),
                    'CHEC': (2, 4),
                    'SCTCam': (1.5, 3)
                    }


@click.command()
@click.argument(
    'input_files', type=click.Path(
        exists=True,
        dir_okay=False,
    ), nargs=-1)
@click.argument(
    'output_file', type=click.Path(
        dir_okay=False,
    ))
@click.option('-n', '--n_events', default=-1)
def main(input_files, output_file, n_events):
    '''
    process multiple simtel files into one hdf containing three groups.
    'runs', 'array_events', 'telescope_events'
    '''

    if os.path.exists(output_file):
        click.confirm(f'File {output_file} exists. Overwrite?', default=False, abort=True)
        os.remove(output_file)

    for input_file in input_files:
        print(f'processing file {input_file}')
        process_file(input_file, output_file, n_events)

    verify_file(output_file)


def process_file(input_file, output_file, n_events=-1):
    event_source = EventSourceFactory.produce(
        input_url=input_file,
        max_events=n_events if n_events > 1 else None,
    )

    calibrator = CameraCalibrator(
        event_source=event_source,
        r1_product='HESSIOR1Calibrator',
    )


    image_features = []
    array_event_information = []
    for event in event_source:
        if number_of_valid_triggerd_cameras(event) >= 2:
            f = calculate_image_features(event, calibrator)
            if len(f) > 1:  # check whtehr at least two telescopes returned hillas features
                array_event_information.append(event_information(event))
                image_features.extend(f)

    df_features = pd.DataFrame(image_features)
    df_features.set_index('telescope_event_id', drop=True, verify_integrity=True, inplace=True)
    fact.io.write_data(df_features, output_file, key='telescope_events', mode='a')

    df_array = pd.DataFrame(array_event_information)
    df_array.set_index('array_event_id', drop=True, verify_integrity=True, inplace=True)
    fact.io.write_data(df_array, output_file, key='array_events', mode='a')

    run_information = read_simtel_mc_information(input_file)

    df_runs = pd.DataFrame([run_information])
    df_runs.set_index('run_id', drop=True, verify_integrity=True, inplace=True)
    fact.io.write_data(df_runs, output_file, key='runs', mode='a')



def verify_file(input_file_path):
    runs = fact.io.read_data(input_file_path, key='runs')
    runs.set_index('run_id', drop=True, verify_integrity=True, inplace=True)

    telescope_events = fact.io.read_data(input_file_path, key='telescope_events')
    telescope_events.set_index('telescope_event_id', drop=True, verify_integrity=True, inplace=True)

    array_events = fact.io.read_data(input_file_path, key='array_events')
    array_events.set_index('array_event_id', drop=True, verify_integrity=True, inplace=True)

    assert len(array_events) == len(telescope_events.array_event_id.unique())
    assert len(runs) == len(telescope_events.run_id.unique())
    assert len(runs) == len(array_events.run_id.unique())

    print(f'Processed {len(runs)} runs, {len(telescope_events)} single telescope events and {len(array_events)} array events.')


def read_simtel_mc_information(simtel_file):
    with pyhessio.open_hessio(simtel_file) as f:
        # do some weird hessio fuckup
        eventstream = f.move_to_next_event()
        _ = next(eventstream)

        d = {
            'mc_spectral_index': f.get_spectral_index(),
            'mc_num_reuse': f.get_mc_num_use(),
            'mc_num_showers': f.get_mc_num_showers(),
            'mc_max_energy': f.get_mc_E_range_Max(),
            'mc_min_energy': f.get_mc_E_range_Min(),
            'mc_max_scatter_range': f.get_mc_core_range_Y(),  # range_X is always 0 in simtel files
            'mc_max_viewcone_radius': f.get_mc_viewcone_Max(),
            'mc_min_viewcone_radius': f.get_mc_viewcone_Min(),
            'run_id': f.get_run_number(),
            'mc_max_altitude': f.get_mc_alt_range_Max(),
            'mc_max_azimuth': f.get_mc_az_range_Max(),
            'mc_min_altitude': f.get_mc_alt_range_Min(),
            'mc_min_azimuth': f.get_mc_az_range_Min(),
        }

        return d


def event_information(event):
    d = {
        'mc_alt': event.mc.alt,
        'mc_az': event.mc.az,
        'mc_core_x': event.mc.core_x,
        'mc_core_y': event.mc.core_y,
        'num_triggered_telescopes': number_of_valid_triggerd_cameras(event),
        'mc_height_first_interaction': event.mc.h_first_int,
        'mc_energy': event.mc.energy.to('TeV').value,
        'mc_corsika_primary_id': event.mc.shower_primary_id,
        'run_id': event.r0.obs_id,
        'array_event_id': generate_unique_array_event_id(event),
    }

    return {k: strip_unit(v) for k, v in d.items()}


def calculate_image_features(event, calibrator):
    '''
    Processes
    '''
    event_info = []
    calibrator.calibrate(event)

    for telescope_id, dl1 in event.dl1.tel.items():
        camera = event.inst.subarray.tels[telescope_id].camera
        if camera.cam_id not in allowed_cameras:
            continue

        telescope_type_name = event.inst.subarray.tels[telescope_id].optics.tel_type
        picture_thresh, boundary_thresh = cleaning_level[camera.cam_id]
        mask = tailcuts_clean(camera, dl1.image[0], boundary_thresh=boundary_thresh, picture_thresh=picture_thresh)

        if mask.sum() < 3:  # only two pixel remaining. No luck anyways.
            continue

        hillas_params = hillas_parameters(
            camera[mask],
            dl1.image[0, mask],
            container=True
        )
        if np.isnan(hillas_params.width.value) or np.isnan(hillas_params.length.value):
            continue


        d = {
            'array_event_id': generate_unique_array_event_id(event),
            'telescope_event_id': generate_unique_telescope_event_id(event, int(telescope_id)),
            'telescope_id': int(telescope_id),
            'camera_name': camera.cam_id,
            'camera_id': names_to_id[camera.cam_id],
            'run_id': event.r0.obs_id,
            'telescope_type_name': telescope_type_name,
            'telescope_type_id': types_to_id[telescope_type_name],
            'pointing_azimuth': event.mc.tel[telescope_id].azimuth_raw,
            'pointing_altitude': event.mc.tel[telescope_id].altitude_raw,
        }

        d.update(hillas_params.as_dict())
        event_info.append({k: strip_unit(v) for k, v in d.items()})

    return event_info


def pair(a, b):
    ''' implmentation of the "elegant pairing function" http://szudzik.com/ElegantPairing.pdf '''
    p = b**2 + a if b > a else a**2 + a + b
    try:
        is_valid_numpy_int = p.dtype in [np.int16, np.int32, np.int64]
        if is_valid_numpy_int:
            return p
    except AttributeError:
        if int(p).bit_length > 63:
            raise ValueError(f'"{a}" and "{b}" cannot be paired. Value exceeds 64 bits.')
        return p


def generate_unique_array_event_id(event):
    return pair(event.r0.obs_id, event.r0.event_id)


def generate_unique_telescope_event_id(event, telescope_id):
    return pair(generate_unique_array_event_id(event), telescope_id)


def number_of_valid_triggerd_cameras(event):
    triggerd_tel_ids = event.trig.tels_with_trigger
    triggerd_camera_names = [event.inst.subarray.tels[i].camera.cam_id for i in triggerd_tel_ids]
    valid_triggered_cameras = list(filter(lambda c: c in allowed_cameras, triggerd_camera_names))
    return len(valid_triggered_cameras)


def strip_unit(v):
    try:
        return v.si.value
    except AttributeError:
        return v


if __name__ == '__main__':
    main()
