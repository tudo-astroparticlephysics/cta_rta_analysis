import click
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
import fact.io


@click.command()
@click.argument('input_dl3_file', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
@click.option('-b', '--bin_edges', nargs=2, help='the lowes and highest bin edges', type=float, default=(0, 0.2))
def main(input_dl3_file, output, bin_edges):
    df = fact.io.read_data(input_dl3_file, key='array_events')

    alt = Angle(df.alt_prediction.values * u.rad).degree
    mc_alt = Angle(df.mc_alt.values * u.rad).degree

    az = Angle(df.az_prediction.values * u.rad).wrap_at(180 * u.deg).degree
    mc_az = Angle(df.mc_az.values * u.rad).wrap_at(180 * u.deg).degree

    distance = np.sqrt((alt - mc_alt)**2 + (az - mc_az)**2)
    # resolution = np.percentile(distance, 68)
    print(f'Plotting a total {len(df)} events')

    plt.hist(distance, bins=np.linspace(*bin_edges, 100), color='blue')
    plt.title('Gammas')
    # plt.axvline(resolution, color='gray', linestyle='--', label='0.68 percentile')
    plt.xlabel('Distance between true and reco position')
    plt.legend()

    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == "__main__":
    main()
