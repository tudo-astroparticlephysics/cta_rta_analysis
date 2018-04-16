import numpy as np
import click
import matplotlib.pyplot as plt
from tqdm import tqdm
import astropy.units as u
from spectrum import CosmicRaySpectrum, CrabSpectrum, MCSpectrum, make_energy_bins, relative_sensitivity
from coordinates import calculate_distance_theta
import fact.io
import pandas as pd
from scipy import optimize
from joblib import Parallel, delayed
import multiprocessing



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


def calculate_sensitivity(
        gammas,
        protons,
        gamma_prediction_mean=0.5,
        signal_region=0.01,
        target_spectrum=CrabSpectrum()
):
    selected_gammas = gammas.query(f'gamma_prediction_mean >={gamma_prediction_mean}').copy()
    selected_protons = protons.query(f'gamma_prediction_mean >={gamma_prediction_mean}').copy()

    if len(selected_gammas) < 10 or len(selected_protons) < 10:
        return np.inf/(u.erg * u.s * u.cm**2)

    n_on, n_off = get_on_and_off_counts(
        selected_protons,
        selected_gammas,
        signal_region_radius=signal_region
    )

    relative_flux = relative_sensitivity(
        n_on,
        n_off,
        alpha=1,
    )

    min_energy = min(gammas.mc_energy.min(),
                     protons.mc_energy.min())
    max_energy = max(gammas.mc_energy.max(),
                     protons.mc_energy.max())

    bin_center = np.mean([min_energy, max_energy]) * u.TeV

    sens = target_spectrum.flux(bin_center) * relative_flux
    return sens


def get_on_and_off_counts(selected_protons, selected_gammas, signal_region_radius):
    """ Get on and off counts from the signal region using a simepl theta**2 cut"""

    # estimate n_off by assuming that the background rate is constant within a
    # smallish theta area around 0. take the mean of the thata square histogram
    # to get a more stable estimate for n_off
    background_region_radius = signal_region_radius
    # # print(selected_protons.theta.count(), selected_protons.weight.count())
    # if selected_protons.weight.count() == 30:
    #     print(selected_protons.weight)
    b = np.nanpercentile(selected_protons.theta**2, 68)

    if b == np.nan:
        b = 10


    H, _ = np.histogram(
        selected_protons.theta**2,
        bins=np.arange(0, b, background_region_radius),
        weights=selected_protons.weight
    )
    n_off = H.mean()

    n_on = selected_gammas.query(f'theta**2 < {signal_region_radius}')['weight'].sum()

    return n_on, n_off


def create_sensitivity_matrix(
            protons,
            gammas,
            n_bins,
            target_spectrum=CrabSpectrum(),
            num_threads=-1
        ):

    min_energy = min(gammas.mc_energy.min(),
                     protons.mc_energy.min())
    max_energy = max(gammas.mc_energy.max(),
                     protons.mc_energy.max())

    edges = np.logspace(np.log10(min_energy), np.log10(
        max_energy), num=n_bins + 1, base=10.0) * u.TeV

    gammas['energy_bin'] = pd.cut(gammas.mc_energy, edges)
    protons['energy_bin'] = pd.cut(protons.mc_energy, edges)

    if num_threads == -1:
        num_threads = multiprocessing.cpu_count()

    if num_threads > 1:
        d = (
            delayed(find_best_sensitivity_in_bin)
            (gammas[gammas.energy_bin == b], protons[protons.energy_bin == b]) for b in gammas.energy_bin.cat.categories
        )

        sensitivity = Parallel(n_jobs=num_threads, verbose=10)(d)
    else:
        sensitivity = [find_best_sensitivity_in_bin(gammas[gammas.energy_bin == b], protons[protons.energy_bin == b]) for b in tqdm(gammas.energy_bin.cat.categories)]



    # multiply the whole thing by the proper unit. There must be a nicer way to do this.
    sensitivity = np.array([s.value for s in sensitivity]) * sensitivity[0].unit
    return sensitivity, edges


def find_best_sensitivity_in_bin(g, p):

    def f(x):
        return calculate_sensitivity(g, p, gamma_prediction_mean=x[0], signal_region=x[1]).value

    ranges = (slice(0.0, 1, 0.025), slice(0.001, 0.08, 0.001))
    # Note: while it seems obviuous to use finish=optimize.fmin here. apparently it
    # tests invalid values. and then everything breaks. Negative theta cuts for
    # example
    res = optimize.brute(f, ranges, finish=None, full_output=True)

    cuts = res[0]
    return calculate_sensitivity(g, p, gamma_prediction_mean=cuts[0], signal_region=cuts[1])



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

    sens, edges = create_sensitivity_matrix(protons, gammas, n_bins, num_threads=n_jobs)
    sens = sens.to(1 / (u.erg * u.s * u.cm**2))
    ax = plot_sensitivity(edges, sens, t_obs, ax=None)
    plot_spectrum(crab, e_min, e_max, ax=ax, color='gray')

    if output:
        plt.savefig(output)
    else:
        plt.show()



if __name__ == '__main__':
    main()
