import click
import matplotlib.pyplot as plt
import fact.io
import astropy.units as u
from spectrum import make_energy_bins


@click.command()
@click.argument('predicted_gammas', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
@click.option('-b', '--bins', default=10, show_default=True)
def main(predicted_gammas, output, bins):
    gammas = fact.io.read_data(predicted_gammas, key='telescope_events')
    gammas_energy_prediction = gammas.groupby('array_event_id')['gamma_energy_prediction'].mean()
    gammas_energy_prediction = gammas_energy_prediction.values * u.TeV

    array_events = fact.io.read_data(predicted_gammas, key='array_events')
    gammas_energy = array_events.mc_energy.values * u.TeV

    bin_edges, _, _ = make_energy_bins(gammas_energy)
    plt.hist(gammas_energy_prediction, bins=bin_edges, histtype='step', lw=2, label='predicted energy')
    plt.hist(gammas_energy, bins=bin_edges, histtype='step', lw=2, label='true energy')
    plt.legend()

    plt.xscale('log')
    plt.xlabel(r'$Energy /  \mathrm{TeV}$')
    plt.ylabel('Counts')

    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == "__main__":
    main()
