import click
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
from astropy.coordinates.angle_utilities import angular_separation
import fact.io


@click.command()
@click.argument('input_files', type=click.Path(exists=True), nargs=-1)
@click.option('-l', '--label', multiple=True)
@click.option('-c', '--color', multiple=True)
@click.option('-o', '--output', type=click.Path(exists=False))
@click.option('-d', '--density', is_flag=True, default=False)
@click.option('-t', '--title', default=None)
def main(input_files, label, color, output, density, title):
    if label and len(input_files) != len(label):
        print('Must pass as many labels as gamma files as proton files')

    if color and label and len(label) != len(color):
        print('Must pass as many colors as labels')

    if not label:
        label = input_files

    if not color:
        color = [None] * len(input_files)

    for input, l, c in zip(input_files, label, color):
        df = fact.io.read_data(input, key='array_events').dropna()

        mc_az = Angle(df.mc_az.values, unit=u.rad).wrap_at(180*u.deg)
        mc_alt = Angle(df.mc_alt.values, unit=u.rad)

        az = Angle(df.az_prediction.values, unit=u.rad).wrap_at(180*u.deg)
        alt = Angle(df.alt_prediction.values, unit=u.rad)

        distance = angular_separation(mc_az, mc_alt, az, alt).to('deg')
        resolution = np.nanpercentile(distance, 68)
        print(f'Plotting a total {len(df)} events')

        plt.hist(distance, bins=np.linspace(0, np.nanpercentile(distance, 95), 100), label=l, color=c, histtype='step', lw=2.5, density=density)

        plt.axvline(resolution, color='gray', linestyle='--', label='0.68 percentile')

    plt.xlabel('Distance between true and reco position')
    if title:
        plt.title(title)
    plt.legend()

    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == "__main__":
    main()
