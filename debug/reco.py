import click
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
import fact.io


def calculate_distance(df):
    alt = Angle(df.alt_prediction.values * u.rad).degree
    mc_alt = Angle(df.mc_alt.values * u.rad).degree

    az = Angle(df.az_prediction.values * u.rad).wrap_at(180 * u.deg).degree
    mc_az = Angle(df.mc_az.values * u.rad).wrap_at(180 * u.deg).degree

    return np.sqrt((alt - mc_alt)**2 + (az - mc_az)**2)


@click.command()
@click.argument('gammas_dl3', type=click.Path(exists=True))
@click.argument('protons_dl3', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
def main(gammas_dl3, protons_dl3, output):

    fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2, 2)

    df_gammas = fact.io.read_data(gammas_dl3, key='array_events')
    distance = calculate_distance(df_gammas)
    ax1.hist(distance, bins=np.linspace(0, 0.2, 50), color='blue')
    ax2.hist(distance**2, bins=np.linspace(0, 0.02, 50), color='blue')


    df_protons = fact.io.read_data(protons_dl3, key='array_events')
    distance = calculate_distance(df_protons)
    ax3.hist(distance, bins=np.linspace(0, 20, 50), color='red')
    ax4.hist(distance**2, bins=np.linspace(0, 40, 50), color='red')

    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == "__main__":
    main()
