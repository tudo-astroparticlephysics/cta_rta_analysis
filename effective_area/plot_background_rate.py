import click
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import fact.io
from fact.analysis import li_ma_significance
from sensitivity import relative_sensitivity
import coordinates
from spectrum import MCSpectrum, CrabSpectrum, CosmicRaySpectrum, make_energy_bins
import pandas as pd


@click.command()
@click.argument('protons_dl3', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
@click.option('--fixed/--no-fixed', default=False)
def main(protons_dl3, output, fixed):
    t_obs = 1 * u.s

    protons = fact.io.read_data(protons_dl3, key='array_events')
    proton_runs = fact.io.read_data(protons_dl3, key='runs')
    mc_production_proton = MCSpectrum.from_cta_runs(proton_runs)

    cosmic = CosmicRaySpectrum()

    protons['weight'] = mc_production_proton.reweigh_to_other_spectrum(cosmic, protons.mc_energy.values * u.TeV, t_assumed_obs=t_obs)

    bin_edges, bin_centers, bin_widths = make_energy_bins(e_min=0.01 * u.TeV, e_max=100 * u.TeV, bins=20)
    protons['energy_bin'] = pd.cut(protons.mc_energy, bin_edges)

    bg = []
    if fixed:
        for bin, p in protons.groupby('energy_bin'):
            print(bin, len(p))
            n = p[p.gamma_prediction_mean > 0.85].weight.sum()
            bg.append(n)

    else:
        for bin, p in protons.groupby('energy_bin'):
            print(bin, len(p))
            _, _, bkg = optimize_background(p)
            bg.append(bkg)

    import IPython; IPython.embed()

    bg = np.array(bg) / bin_widths.to(u.MeV) / ((1 - np.cos(12 * u.deg)) * u.sr) /u.s
    plt.step(bin_centers, bg, where='mid')
    plt.xscale('log')
    plt.yscale('log')
    if output:
        plt.savefig(output)
    else:
        plt.show()


def optimize_background(protons, min_events=10):
    if len(protons) < 10:
        return np.nan, np.nan

    cuts = np.linspace(0, 1, 100)

    min_cut = 0
    min_bkg = np.inf
    for c in cuts:
        # print(c, len(protons[protons.gamma_prediction_mean > c]))
        if len(protons[protons.gamma_prediction_mean > c]) < 10:
            continue

        n = protons[protons.gamma_prediction_mean > c].weight.sum()
        if n < min_bkg:
            min_bkg = n
            min_cut = c

    return min_cut, len(protons[protons.gamma_prediction_mean > min_cut]), min_bkg


if __name__ == '__main__':
    main()
