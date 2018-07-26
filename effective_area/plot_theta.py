import click
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import fact.io
from fact.analysis import li_ma_significance
from sensitivity import relative_sensitivity
import coordinates
from spectrum import MCSpectrum, CrabSpectrum, CosmicRaySpectrum


@click.command()
@click.argument('gammas_dl3', type=click.Path(exists=True))
@click.argument('protons_dl3', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
def main(gammas_dl3, protons_dl3, output):

    t_obs = 0.5 * u.h

    gammas = fact.io.read_data(gammas_dl3, key='array_events')

    gamma_runs = fact.io.read_data(gammas_dl3, key='runs')
    mc_production_gamma = MCSpectrum.from_cta_runs(gamma_runs)

    protons = fact.io.read_data(protons_dl3, key='array_events')

    e_min = gammas.mc_energy.min() * u.TeV
    e_max = gammas.mc_energy.max() * u.TeV

    gammas['theta'] = coordinates.calculate_distance_theta(gammas, source_az=0 * u.deg, source_alt=70 * u.deg).to(u.deg).value
    protons['theta'] = coordinates.calculate_distance_theta(protons, source_az=0 * u.deg, source_alt=70 * u.deg).to(u.deg).value

    # gammas = gammas[gammas.num_triggered_telescopes >= 4]
    # protons = protons[protons.num_triggered_telescopes >= 4]

    proton_runs = fact.io.read_data(protons_dl3, key='runs')
    mc_production_proton = MCSpectrum.from_cta_runs(proton_runs)

    crab = CrabSpectrum()
    cosmic = CosmicRaySpectrum()

    gammas['weight'] = mc_production_gamma.reweigh_to_other_spectrum(crab, gammas.mc_energy.values * u.TeV, t_assumed_obs=t_obs)
    protons['weight'] = mc_production_proton.reweigh_to_other_spectrum(cosmic, protons.mc_energy.values * u.TeV, t_assumed_obs=t_obs)

    # gammas['weight'] = 1
    # protons['weight'] = 1

    significance = 0
    best_prediction_cut = 0
    best_theta_cut = 0
    for cut in np.linspace(0., 1, 15):
        for theta_cut in np.linspace(0.01, 0.3, 10):
            gammas_gammalike = gammas.query(f'gamma_prediction_mean > {cut}').copy()
            protons_gammalike = protons.query(f'gamma_prediction_mean > {cut}').copy()

            on = gammas_gammalike
            off = protons_gammalike

            off_bins = np.arange(0, 1, theta_cut**2)
            h, _ = np.histogram(off['theta']**2, bins=off_bins, weights=off.weight)

            n_off = h.mean()
            n_on = on.query(f'theta <= {theta_cut}').weight.sum() + n_off


            s = li_ma_significance(n_on, n_off, alpha=1)
            print(f'significance: {s}, prediction cut: {cut}, theta cut: {theta_cut}')
            if s > significance:
                significance = s
                best_prediction_cut = cut
                best_theta_cut = theta_cut

    # best_prediction_cut = 0.9
    gammas_gammalike = gammas.query(f'gamma_prediction_mean > {best_prediction_cut}').copy()
    protons_gammalike = protons.query(f'gamma_prediction_mean > {best_prediction_cut}').copy()

    print(f'Cut: {best_prediction_cut} {best_theta_cut}, gammas: {len(gammas_gammalike)}, protons: {len(protons_gammalike)}')

    on = gammas_gammalike
    off = protons_gammalike

    bins = np.arange(0, 1, 0.01)
    h_off, _ = np.histogram(off['theta']**2, bins=bins, weights=off.weight)
    h_on, _ = np.histogram(on['theta']**2, bins=bins, weights=on.weight)

    plt.step(bins[:-1], h_on + h_off.mean(), where='post', label='on events')
    plt.step(bins[:-1], h_off, where='post', label='off events')
    plt.ylim([0, max(h_on + h_off) * 1.3])

    n_off = h_off.mean()
    n_on = on.query(f'theta <= {best_theta_cut}').weight.sum() + n_off

    sens = relative_sensitivity(n_on, n_off, alpha=1)
    rate = crab.rate(e_min, e_max)
    print(f'Integral sensitivity {rate.to("1/(cm^2 s)") * sens}, relative flux {sens}')

    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == '__main__':
    main()
