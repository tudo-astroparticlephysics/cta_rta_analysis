import click
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import Angle
import fact.io
import numpy as np
from matplotlib.colors import LogNorm
# from coordinates import horizontal_to_skycoord, wrap_angles


@click.command()
@click.argument('gammas_dl3', type=click.Path(exists=True))
@click.argument('protons_dl3', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
def main(gammas_dl3, protons_dl3, output):


    gammas = fact.io.read_data(gammas_dl3, key='array_events')
    protons = fact.io.read_data(protons_dl3, key='array_events')

    fig, ([ax1, ax2], [ax3, ax4]) = plt.subplots(2, 2, figsize=(10, 10), constrained_layout=True)

    az = Angle(protons.mc_az.values, unit=u.rad).wrap_at(180 * u.deg).degree
    alt = Angle(protons.mc_alt.values, unit=u.rad).degree

    # ax1.scatter(ra.degree, dec.degree, label='proton', s=2, alpha=0.5)
    ax1.scatter(alt, az, label=f'proton ({len(az)})', s=2, alpha=0.5)

    az = Angle(gammas.mc_az.values, unit=u.rad).wrap_at(180 * u.deg).degree
    alt = Angle(gammas.mc_alt.values, unit=u.rad).degree

    # ax1.scatter(ra.degree, dec.degree, label='gamma', s=2, alpha=0.5)
    ax1.scatter(alt, az, label=f'gamma ({len(az)})', s=8, alpha=0.5)
    ax1.set_title('True Direction')



    az = Angle(protons.az_prediction.values, unit=u.rad).wrap_at(180 * u.deg).degree
    alt = Angle(protons.alt_prediction.values, unit=u.rad).degree

    ax2.scatter(alt, az, label=f'proton ({len(az)})', s=2, alpha=0.5)


    az = Angle(gammas.az_prediction.values, unit=u.rad).wrap_at(180 * u.deg).degree
    alt = Angle(gammas.alt_prediction.values, unit=u.rad).degree

    ax2.scatter(alt, az, label=f'gamma ({len(az)})', s=2, alpha=0.5)
    ax2.set_title('Reconstruncted Direction')

    # ax1.set_ylim([-25, 25])
    # ax1.set_xlim([45, 95])
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_ylim(ax1.get_ylim())

    ax2.set_xlabel('alt')
    ax2.set_ylabel('az')

    ax1.set_xlabel('alt')
    ax1.set_ylabel('az')


    az = Angle(gammas.az_prediction.values, unit=u.rad).wrap_at(180 * u.deg).degree
    alt = Angle(gammas.alt_prediction.values, unit=u.rad).degree

    az = np.nan_to_num(az)
    alt = np.nan_to_num(alt)

    bins = [np.linspace(*ax1.get_xlim(), 100), np.linspace(*ax1.get_ylim(), 100)]
    ax3.hist2d(alt, az, bins=bins, label='gammas_prediction', norm=LogNorm())
    ax3.set_ylim(ax1.get_ylim())
    ax3.set_xlim(ax1.get_xlim())

    az = Angle(protons.az_prediction.values, unit=u.rad).wrap_at(180 * u.deg).degree
    alt = Angle(protons.alt_prediction.values, unit=u.rad).degree

    az = np.nan_to_num(az)
    alt = np.nan_to_num(alt)
    # import IPython; IPython.embed()

    ax4.hist2d(alt, az, bins=bins, label='proton_prediction', norm=LogNorm())
    # ax4.set_ylim(ax1.get_ylim())
    # ax4.set_xlim(ax1.get_xlim())

    ax1.legend()
    ax2.legend()

    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == '__main__':
    main()
