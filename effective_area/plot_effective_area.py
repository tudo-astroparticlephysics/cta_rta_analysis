import click
import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import binom_conf_interval
import fact.io
from spectrum import MCSpectrum, make_energy_bins
import astropy.units as u


@click.command()
@click.argument('gammas_dl3', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
@click.option('-b', '--bins', default=10, show_default=True)
def main(gammas_dl3, output, bins):
    gammas = fact.io.read_data(gammas_dl3, key='array_events')
    gammas_energy = gammas.mc_energy.values

    runs = fact.io.read_data(gammas_dl3, key='runs')
    mc_production = MCSpectrum.from_cta_runs(runs)
    bins, bin_center, bin_widths = make_energy_bins(gammas.mc_energy.values * u.TeV, bins=40)

    hist_all = mc_production.expected_events_for_bins(energy_bins=bins)
    hist_selected, _ = np.histogram(gammas_energy, bins=bins)

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

    plt.errorbar(
        bin_center.value,
        area.value,
        xerr=bin_widths.value / 2.0,
        yerr=[lower.value, upper.value],
        linestyle='',
    )

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$E_{\mathrm{True}} /  \mathrm{TeV}$')
    plt.ylabel(r'$\mathrm{Mean Effective\; Area} / \mathrm{m}^2$')

    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == "__main__":
    main()
