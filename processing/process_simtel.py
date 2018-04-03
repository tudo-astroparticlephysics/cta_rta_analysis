from ctapipe.io.eventsourcefactory import EventSourceFactory
from ctapipe.calib import CameraCalibrator
from ctapipe.image.hillas import hillas_parameters
from ctapipe.image.cleaning import tailcuts_clean
import pandas as pd
import fact.io
import click
from tqdm import tqdm
import numpy as np
import os
import pyhessio


names_to_id = {'LSTCam': 1, 'NectarCam': 2, 'FlashCam': 3, 'DigiCam': 4, 'CHEC': 5}
types_to_id = {'LST': 1, 'MST': 2, 'SST': 3}
allowed_cameras = ['LSTCam', 'NectarCam', 'DigiCam']


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

    run_information = []
    count = 0
    for input_file in input_files:
        print(f'processing file {input_file}')
        run_information.append(read_simtel_mc_information(input_file))
        count += process_file(input_file, output_file, n_events)

    df_runs = pd.DataFrame(run_information)
    df_runs.set_index('run_id', drop=True, verify_integrity=True, inplace=True)
    fact.io.write_data(df_runs, output_file, key='runs')

    print(f'Processed {count} (telescope-wise) events in total.')


def process_file(input_file, output_file, n_events):
    event_source = EventSourceFactory.produce(
        input_url=input_file,
        max_events=n_events if n_events > 1 else None,
    )

    calibrator = CameraCalibrator(
        event_source=event_source,
    )


    image_features = []
    array_event_information = []
    for event in tqdm(event_source):
        if number_of_valid_triggerd_cameras(event) >= 2:
            array_event_information.append(event_information(event))
            image_features.extend(calculate_image_features(event, calibrator))

    df_features = pd.DataFrame(image_features)
    df_features.set_index('telescope_event_id', drop=True, verify_integrity=True, inplace=True)
    fact.io.write_data(df_features, output_file, key='telescope_events', mode='a')

    df_array = pd.DataFrame(array_event_information)
    df_array.set_index('array_event_id', drop=True, verify_integrity=True, inplace=True)
    fact.io.write_data(df_array, output_file, key='array_events', mode='a')
    return len(df_features)


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

        mask = tailcuts_clean(camera, dl1.image[0], min_number_picture_neighbors=2)
        hillas_params = hillas_parameters(
            camera[mask],
            dl1.image[0, mask],
            container=True
        )

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
            'pointing_altitude': ((np.pi / 2) - event.mc.tel[telescope_id].altitude_raw),
        }

        d.update(hillas_params.as_dict())
        event_info.append({k: strip_unit(v) for k, v in d.items()})

    return event_info


def generate_unique_array_event_id(event):
    return (event.r0.obs_id << 22) + event.r0.event_id + event.count


def generate_unique_telescope_event_id(event, telescope_id):
    return (event.r0.obs_id << 10) + (event.r0.event_id << 8) + (event.count << 4) + telescope_id


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
