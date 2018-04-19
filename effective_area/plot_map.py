import click
import matplotlib.pyplot as plt
import astropy.units as u
import fact.io
# from coordinates import horizontal_to_skycoord, wrap_angles


@click.command()
@click.argument('gammas_dl3', type=click.Path(exists=True))
@click.argument('protons_dl3', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
def main(gammas_dl3, protons_dl3, output):


    gammas = fact.io.read_data(gammas_dl3, key='array_events')
    protons = fact.io.read_data(protons_dl3, key='array_events')


    fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(10, 5))

    az = protons.mc_az.values * u.rad
    alt = protons.mc_alt.values * u.rad

    # ax1.scatter(ra.degree, dec.degree, label='proton', s=2, alpha=0.5)
    ax1.scatter(*wrap_angles(alt, az), label=f'proton ({len(az)})', s=2, alpha=0.5)

    az = gammas.mc_az.values * u.rad
    alt = gammas.mc_alt.values * u.rad

    # ax1.scatter(ra.degree, dec.degree, label='gamma', s=2, alpha=0.5)
    ax1.scatter(*wrap_angles(alt, az), label=f'gamma ({len(az)})', s=8, alpha=0.5)
    ax1.set_title('True Direction')


    az = protons.az_prediction.values * u.rad
    alt = protons.alt_prediction.values * u.rad

    ax2.scatter(*wrap_angles(alt, az), label=f'proton ({len(az)})', s=2, alpha=0.5)

    az = gammas.az_prediction.values * u.rad
    alt = gammas.alt_prediction.values * u.rad

    ax2.scatter(*wrap_angles(alt, az), label=f'gamma ({len(az)})', s=2, alpha=0.5)
    ax2.set_title('Reconstruncted Direction')

    ax1.set_ylim([-25, 25])
    ax1.set_xlim([45, 95])
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_ylim(ax1.get_ylim())

    ax2.set_xlabel('alt')
    ax2.set_ylabel('az')

    ax1.set_xlabel('alt')
    ax1.set_ylabel('az')


    ax1.legend()
    ax2.legend()

    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == '__main__':
    main()
