import click
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
import fact.io
import coordinates
from spectrum import MCSpectrum, CrabSpectrum, CosmicRaySpectrum, make_energy_bins
from fact.analysis import li_ma_significance


@click.command()
@click.argument('gammas', type=click.Path(exists=True))
@click.argument('protons', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
def main(gammas, protons, output):

    t_obs = 50 * u.h

    gammas = fact.io.read_data(gammas, key='array_events')
    gammas = gammas.dropna()


    gamma_runs = fact.io.read_data(gammas, key='runs')
    mc_production_gamma = MCSpectrum.from_cta_runs(gamma_runs)

    protons = fact.io.read_data(protons, key='array_events')
    protons = protons.dropna()

    # print(f'Plotting {len(protons)} protons and {len(gammas)} gammas.')
    proton_runs = fact.io.read_data(protons, key='runs')
    mc_production_proton = MCSpectrum.from_cta_runs(proton_runs)


    crab = CrabSpectrum()
    cosmic = CosmicRaySpectrum()


    gammas['weight'] = mc_production_gamma.reweigh_to_other_spectrum(crab, gammas.mc_energy.values * u.TeV, t_assumed_obs=t_obs)
    protons['weight'] = mc_production_proton.reweigh_to_other_spectrum(cosmic, protons.mc_energy.values * u.TeV, t_assumed_obs=t_obs)


    # gammas_gammalike = gammas.query(f'gamma_prediction_mean > {cut}')
    # protons_gammalike = protons.query(f'gamma_prediction_mean > {cut}')


    bin_edges, _, _ = make_energy_bins(gammas.mc_energy.values * u.TeV, bins=20)
    on, off, alpha = coordinates.split_on_off(gammas, protons, on_region_radius=0.4 * u.deg)
    print(f'alpha:{alpha}')
    on['energy_bin'] = pd.cut(on.mc_energy, bin_edges)
    off['energy_bin'] = pd.cut(off.mc_energy, bin_edges)
    for ((_, g_on), (_, g_off)) in zip(on.groupby('energy_bin'), off.groupby('energy_bin')):
        n_on = g_on.weight.sum()
        n_off = g_off.weight.sum()
        print('----'*20)
        print(n_on, n_off)
        print(g_on.size, g_off.size)
        print(li_ma_significance(n_on, n_off, alpha=1))

    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == '__main__':
    main()
