import click
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
from astropy.coordinates.angle_utilities import angular_separation
import fact.io


@click.command()
@click.argument('input_dl3_file', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
@click.option('-t', '--title', default=None)
def main(input_dl3_file, output, title):
    df = fact.io.read_data(input_dl3_file, key='array_events').dropna()

    mc_az = Angle(df.mc_az.values, unit=u.rad).wrap_at(180*u.deg).degree
    mc_alt = Angle(df.mc_alt.values, unit=u.rad).degree

    az = Angle(df.az_prediction.values, unit=u.rad).degree
    alt = Angle(df.alt_prediction.values, unit=u.rad).degree

    fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(7, 10))
    bins = [np.linspace(60, 80, 60), np.linspace(-10, 10, 60)]
    ax1.hist2d(mc_alt, mc_az, bins=bins, cmap='gray')
    ax1.set_xlabel('True Alt')
    ax1.set_ylabel('True Az')
    ax1.set_aspect('equal')

    ax2.hist2d(alt, az, bins=bins, cmap='gray')
    ax2.set_aspect('equal')
    ax2.set_xlabel('Predicted Alt')
    ax2.set_ylabel('Predicted Az')

    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == "__main__":
    main()
