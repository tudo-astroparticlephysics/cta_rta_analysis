import click
import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import binom_conf_interval
import fact.io
from spectrum import MCSpectrum, make_energy_bins
import astropy.units as u


@click.command()
@click.option('-g', '--gamma_input', type=click.Path(exists=True), multiple=True)
@click.option('-l', '--label', multiple=True)
@click.option('-c', '--color', multiple=True)
@click.option('-o', '--output', type=click.Path(exists=False))
@click.option('-b', '--n_bins', default=20, show_default=True)
@click.option('-t', '--threshold', default=0.0, show_default=True, help='prediction threshold to apply')
def main(gamma_input, label, color, output, n_bins, threshold):

    if label and len(gamma_input) != len(label):
        print('Must pass as many labels as gamma files as proton files')

    if color and len(label) != len(color):
        print('Must pass as many colors as labels')

    if not label:
        label = [None] * len(gamma_input)

    if not color:
        color = [None] * len(gamma_input)

    bins, bin_center, bin_widths = make_energy_bins(e_min=0.003 * u.TeV, e_max=300 * u.TeV, bins=n_bins)

    for input, l, c in zip(gamma_input, label, color):
        print(f'reading file {input} for label {l} with color {c}')

        gammas = fact.io.read_data(input, key='array_events').dropna()

        if threshold > 0:
            gammas = gammas.loc[gammas.gamma_prediction_mean >= threshold]

        gammas_energy = gammas.mc_energy.values

        runs = fact.io.read_data(input, key='runs')
        mc_production = MCSpectrum.from_cta_runs(runs)

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
            color=c,
            label=l
        )


    plt.ylim([10, 10E7])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$E_{\mathrm{True}} /  \mathrm{TeV}$')
    plt.ylabel(r'$\mathrm{Mean Effective\; Area} / \mathrm{m}^2$')

    if label[0]:
        plt.legend()

    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == "__main__":
    main()
