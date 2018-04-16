import numpy as np
import click
import matplotlib.pyplot as plt
import astropy.units as u
from spectrum import CosmicRaySpectrum, CrabSpectrum, MCSpectrum, make_energy_bins
from coordinates import calculate_distance_theta
import fact.io
import pandas as pd
from sensitivity import find_differential_sensitivity


@u.quantity_input(bin_edges=u.TeV, t_obs=u.h)
def plot_sensitivity(bin_edges, sensitivity, t_obs, ax=None, scale=True, **kwargs):
    error = None

    if not ax:
        _, ax = plt.subplots(1)

    bin_center = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width = np.diff(bin_edges)

    sensitivity = sensitivity.to(1 / (u.erg * u.s * u.cm**2))

    if sensitivity.ndim == 2:
        error = sensitivity.std(axis=0) / 2
        sensitivity = sensitivity.mean(axis=0)

    if scale:
        sensitivity = sensitivity * bin_center.to('erg')**2
        if error:
            error = error * bin_center.to('erg')**2

    ax.errorbar(
        bin_center.value,
        sensitivity.value,
        xerr=bin_width.value * 0.5,
        yerr=error.value if error else None,
        marker='.',
        linestyle='',
        capsize=0,
        **kwargs,
    )

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel(r'$ \mathrm{photons} / \mathrm{erg s} \mathrm{cm}^2$ in ' + str(t_obs.to('h')) + ' hours' )
    if scale:
        ax.set_ylabel(r'$ E^2 \cdot \mathrm{photons} \quad \mathrm{erg} /( \mathrm{s} \quad  \mathrm{cm}^2$ )  in ' + str(t_obs.to('h')) )
    ax.set_xlabel(r'$E /  \mathrm{TeV}$')

    return ax


@u.quantity_input(e_min=u.TeV, e_max=u.TeV)
def plot_spectrum(spectrum, e_min, e_max, ax=None, scale=True, **kwargs):

    if not ax:
        _, ax = plt.subplots(1)

    e = np.linspace(e_min, e_max, 1000)
    flux = spectrum.flux(e).to(1 / (u.erg * u.s * u.cm**2))

    if scale:
        flux = flux * e.to('erg')**2

    ax.plot(
        e,
        flux,
        linestyle='--',
        **kwargs,
    )
    return ax




@click.command()
@click.argument('gammas_dl3', type=click.Path(exists=True, dir_okay=False,))
@click.argument('protons_dl3', type=click.Path(exists=True, dir_okay=False,))
@click.option('-n', '--n_bins', type=click.INT, default=10, help='energy bins to plot')
@click.option('-j', '--n_jobs', type=click.INT, default=-1, help='number of threads to use inparallel')
@click.option('-o', '--output', type=click.Path(exists=False))
def main(
    gammas_dl3,
    protons_dl3,
    n_bins,
    n_jobs,
    output,
):
    '''
    Plots a sensitivity curve vs real energy. For each energy bin it performs a gridsearch
    to find the theta and gamma_prediction_mean cuts that produce the highest sensitivity.
    '''
    t_obs = 50 * u.h
    e_min, e_max = 0.003 * u.TeV, 300 * u.TeV
    bin_edges, _, _ = make_energy_bins(e_min=e_min, e_max=e_max, bins=n_bins)


    columns = ['gamma_prediction_mean', 'az_prediction', 'alt_prediction', 'mc_alt', 'mc_az', 'mc_energy']
    gammas = fact.io.read_data(gammas_dl3, key='array_events', columns=columns)
    gammas = gammas.dropna()


    gamma_runs = fact.io.read_data(gammas_dl3, key='runs')
    mc_production_gamma = MCSpectrum.from_cta_runs(gamma_runs)

    protons = fact.io.read_data(protons_dl3, key='array_events', columns=columns)
    protons = protons.dropna()

    # print(f'Plotting {len(protons)} protons and {len(gammas)} gammas.')
    proton_runs = fact.io.read_data(protons_dl3, key='runs')
    mc_production_proton = MCSpectrum.from_cta_runs(proton_runs)


    crab = CrabSpectrum()
    cosmic = CosmicRaySpectrum()


    gammas['weight'] = mc_production_gamma.reweigh_to_other_spectrum(crab, gammas.mc_energy.values * u.TeV, t_assumed_obs=t_obs)
    protons['weight'] = mc_production_proton.reweigh_to_other_spectrum(cosmic, protons.mc_energy.values * u.TeV, t_assumed_obs=t_obs)

    protons['energy_bin'] = pd.cut(protons.mc_energy, bin_edges)
    gammas['energy_bin'] = pd.cut(gammas.mc_energy, bin_edges)

    gammas['theta'] = calculate_distance_theta(gammas)
    protons['theta'] = calculate_distance_theta(protons)

    sens, edges = find_differential_sensitivity(protons, gammas, n_bins, num_threads=n_jobs)
    sens = sens.to(1 / (u.erg * u.s * u.cm**2))
    ax = plot_sensitivity(edges, sens, t_obs, ax=None)
    plot_spectrum(crab, e_min, e_max, ax=ax, color='gray')

    if output:
        plt.savefig(output)
    else:
        plt.show()



if __name__ == '__main__':
    main()
