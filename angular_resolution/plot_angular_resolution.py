import click
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
from matplotlib.colors import LogNorm
import fact.io


@click.command()
@click.argument('input_dl3_file', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
def main(input_dl3_file, output):
    df = fact.io.read_data(input_dl3_file, key='events').dropna()

    alt = Angle(df.alt_prediction.values * u.rad).degree
    mc_alt = Angle(df.mc_alt.values * u.rad).degree

    az = Angle(df.az_prediction.values * u.rad).wrap_at(180 * u.deg).degree
    mc_az = Angle(df.mc_az.values * u.rad).wrap_at(180 * u.deg).degree

    distance = np.sqrt((alt - mc_alt)**2 + (az - mc_az)**2)
    resolution = np.percentile(distance, 68)

    fig, (ax1, ax2) = plt.subplots(2, 1)

    (_, _, _, im) = ax1.hist2d(
        alt, az, range=[[69.5, 70.5], [-0.4, 0.4]], bins=200, cmap='viridis',
    )

    ax1.set_ylabel('Azimuth')
    ax1.set_ylim([-0.2, 0.2])
    ax1.set_xlim([69.8, 70.2])
    ax1.get_xaxis().set_visible(False)
    ax1.grid(b=False)
    ax1.plot(np.unique(mc_alt), np.unique(mc_az), '+', ms=20, mew=1, color='red')
    fig.colorbar(im, ax=ax1)

    (_, _, _, im) = ax2.hist2d(
        alt, az, range=[[69.5, 70.5], [-0.4, 0.4]], bins=200, cmap='viridis', norm=LogNorm()
    )
    ax2.set_xlabel('Altitude')
    ax2.set_ylim([-0.2, 0.2])
    ax2.set_xlim([69.8, 70.2])
    ax2.grid(b=False)
    ax2.plot(np.unique(mc_alt), np.unique(mc_az), '+', ms=20, mew=1, color='red')
    props = dict(facecolor='white', alpha=0.7)
    ax2.text(69.85,
             -0.17,
             'Angular Resolution ~ {:.3f} degrees'.format(resolution),
             bbox=props
             )
    fig.colorbar(im, ax=ax2)

    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == "__main__":
    main()
