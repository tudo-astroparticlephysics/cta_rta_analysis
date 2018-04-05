import click
import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import binom_conf_interval
import fact.io
from spectrum import MCSpectrum


@click.command()
@click.argument('predicted_gammas', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
@click.option('-b', '--bins', default=10, show_default=True)
def main(predicted_gammas, output, bins):
    gammas = fact.io.read_data(predicted_gammas, key='array_events')
    gammas_energy = gammas.mc_energy.values

    runs = fact.io.read_data(predicted_gammas, key='runs')
    mc_production = MCSpectrum.from_cta_runs(runs)

    hist_all, bin_edges = mc_production.expected_events_for_bins(bins=10)
    hist_selected, _ = np.histogram(gammas_energy, bins=bin_edges)

    invalid = hist_selected > hist_all
    hist_selected[invalid] = hist_all[invalid]
    # use astropy to compute errors on that stuff
    lower_conf, upper_conf = binom_conf_interval(hist_selected, hist_all)

    # scale confidences to match and split
    lower_conf = lower_conf * mc_production.generation_area
    upper_conf = upper_conf * mc_production.generation_area

    area = (hist_selected / hist_all) * mc_production.generation_area

    # matplotlib wants relative offsets for errors. the conf values are absolute.
    lower = area - lower_conf
    upper = upper_conf - area


    bin_center = 0.5 * (bin_edges[:-1].value + bin_edges[1:].value)
    bin_width = np.diff(bin_edges).value

    plt.errorbar(
        bin_center,
        area.value,
        xerr=bin_width / 2.0,
        yerr=[lower.value, upper.value],
        linestyle='',
    )

    plt.xscale('log')
    plt.xlabel(r'$E_{\mathrm{True}} /  \mathrm{TeV}$')
    plt.ylabel(r'$\mathrm{Mean Effective\; Area} / \mathrm{m}^2$')

    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == "__main__":
    main()
