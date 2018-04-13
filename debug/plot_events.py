
from ctapipe.io.eventsourcefactory import EventSourceFactory
from ctapipe.calib import CameraCalibrator
from ctapipe.image.hillas import hillas_parameters
from ctapipe.image.cleaning import tailcuts_clean
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib.pyplot as plt
from ctapipe.visualization import CameraDisplay
from tqdm import tqdm
import click


names_to_id = {'LSTCam': 1, 'NectarCam': 2, 'FlashCam': 3, 'DigiCam': 4, 'CHEC': 5}
types_to_id = {'LST': 1, 'MST': 2, 'SST': 3}
allowed_cameras = ['LSTCam', 'NectarCam', 'DigiCam']


cleaning_level = {
                    'ASTRICam': (5, 7),  # (5, 10)?
                    'FlashCam': (12, 15),
                    'LSTCam': (3.5, 6),  # ?? (3, 6) for Abelardo...
                    # ASWG Zeuthen talk by Abelardo Moralejo:
                    'NectarCam': (4, 8),
                    # "FlashCam": (4, 8),  # there is some scaling missing?
                    'DigiCam': (3, 6),
                    'CHEC': (2, 4),
                    'SCTCam': (1.5, 3)
                    }


@click.command()
@click.argument('simtel_file', type=click.Path(exists=True))
@click.argument('output_pdf', type=click.Path(exists=False))
@click.option('-n', '--num_events', default=10)
def main(simtel_file, output_pdf, num_events):

    event_source = EventSourceFactory.produce(
        input_url=simtel_file,
        max_events=num_events,
    )

    calibrator = CameraCalibrator(
        eventsource=event_source,
        # r1_product='HESSIOR1Calibrator',
    )
    with PdfPages(output_pdf) as pdf:
        for event in tqdm(event_source, total=num_events):
            calibrator.calibrate(event)

            plt.figure(figsize=(24, 24))
            plt.suptitle(f'EVENT {event.r0.event_id} \n Energy: {event.mc.energy} \n Type: {event.mc.shower_primary_id}')
            plot_event(event, pdf)
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()


def plot_event(event, pdf):

    cams = [event.inst.subarray.tels[i].camera for i in event.r0.tels_with_data]
    cams = [c for c in cams if c.cam_id in allowed_cameras]
    n_tels = len(cams)

    p = 1
    for telescope_id, dl1 in event.dl1.tel.items():
        camera = event.inst.subarray.tels[telescope_id].camera
        if camera.cam_id not in allowed_cameras:
            continue

        nn = int(np.ceil(np.sqrt(n_tels)))
        ax = plt.subplot(nn, nn, p)
        p +=1

        boundary_thresh, picture_thresh = cleaning_level[camera.cam_id]
        mask = tailcuts_clean(camera, dl1.image[0], boundary_thresh=boundary_thresh, picture_thresh=picture_thresh, min_number_picture_neighbors=1)

        if mask.sum() < 3:  # only two pixel remaining. No luck anyways.
            continue

        hillas_params = hillas_parameters(
            camera[mask],
            dl1.image[0, mask],
            container=False
        )

        disp = CameraDisplay(camera, ax=ax, title="CT{0}".format(telescope_id))
        disp.pixels.set_antialiaseds(False)
        disp.autoupdate = False
        disp.add_colorbar()

        # Show the camera image and overlay Hillas ellipse and clean pixels
        disp.image = dl1.image[0]
        disp.cmap = 'viridis'
        disp.highlight_pixels(mask, color='black')
        disp.overlay_moments(hillas_params, color='cyan', linewidth=3)


if __name__ == '__main__':
    main()
