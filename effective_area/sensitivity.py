from scipy.optimize import minimize_scalar
from fact.analysis import li_ma_significance
import astropy.units as u
import numpy as np
import pandas as pd
from spectrum import CrabSpectrum
import multiprocessing
from scipy import optimize
from tqdm import tqdm
from astropy.table import QTable, Table
import coordinates
from itertools import starmap, zip_longest



def relative_sensitivity(
        n_on,
        n_off,
        alpha,
        target_significance=5,
        significance_function=li_ma_significance,
        ):
    '''
    Calculate the relative sensitivity defined as the flux
    relative to the reference source that is detectable with
    significance in t_ref.

    Given measured `n_on` and `n_off` during a time period `t_obs`,
    we estimate the number of gamma events `n_signal` as `n_on - alpha * n_off`.
    The number of background events `n_background` is estimated as
    `n_off * alpha`.

    So we find the relative sensitivity as the scaling factor for `n_signal`
    that yields a significance of `target_significance`.


    Parameters
    ----------
    n_on: int or array-like
        Number of signal-like events for the on observations
    n_off: int or array-like
        Number of signal-like events for the off observations
    alpha: float
        Scaling factor between on and off observations.
        1 / number of off regions for wobble observations.
    t_obs: astropy.units.Quantity of type time
        Total observation time
    t_ref: astropy.units.Quantity of type time
        Reference time for the detection
    significance: float
        Significance necessary for a detection
    significance_function: function
        A function f(n_on, n_off, alpha) -> significance in sigma
        Used to calculate the significance, default is the Li&Ma
        likelihood ratio test formula.
        Li, T-P., and Y-Q. Ma.
        "Analysis methods for results in gamma-ray astronomy."
        The Astrophysical Journal 272 (1983): 317-324.
        Formula (17)
    '''
    if np.isnan(n_on) or np.isnan(n_off):
        return np.nan

    if n_off * alpha > n_on:
        return np.nan

    n_background = n_off * alpha
    n_signal = n_on - n_background



    def equation(relative_flux):
        n_on_scaled = n_signal * relative_flux + n_background
        s = significance_function(n_on_scaled, n_background, alpha=1) - target_significance
        # print(n_on, n_off, relative_flux, s, significance_function(n_on, n_off, alpha))
        return s
    try:
        phi_rel = optimize.newton(equation, x0=0.001)
    except RuntimeError:
        # warnings.warn('Could not calculate relative significance, returning nan')
        phi_rel = np.nan

    s = significance_function(n_signal * phi_rel + n_background, n_off, alpha)
    if s > 5.1 or s < 4.9:
        print('did the thing effing fail?')
    return phi_rel


relative_sensitivity = np.vectorize(
    relative_sensitivity,
    excluded=[
        't_obs',
        't_ref',
        'alpha',
        'target_significance',
        'significance_function',
    ]
)



# def relative_sensitivity(
#         n_on,
#         n_off,
#         alpha,
#         target_significance=5,
# ):
#     '''
#     Calculate the relative sensitivity defined as the flux
#     relative to the reference source that is detectable with
#     significance in t_ref.
#
#     Parameters
#     ----------
#     n_on: int or array-like
#         Number of signal-like events for the on observations
#     n_off: int or array-like
#         Number of signal-like events for the off observations
#     alpha: float
#         Scaling factor between on and off observations.
#         1 / number of off regions for wobble observations.
#     target_significance: float
#         Significance necessary for a detection
#
#     Returns
#     ----------
#     The relative flux neccessary to detect the source with the given target significance.
#     '''
#     is_scalar = np.isscalar(n_on) and np.isscalar(n_off)
#
#     if is_scalar:
#         n_on = [n_on]
#         n_off = [n_off]
#
#     scale = []
#     for on, off in zip(n_on, n_off):
#         def f(relative_flux):
#
#             s = li_ma_significance((on - off) * relative_flux + off, off, alpha=alpha)
#             print(relative_flux, s)
#             return (target_significance - s)**2
#
#         s = minimize_scalar(f, bounds=(1e-13, 13000), method='bounded')
#
#         scale.append(s.x)
#
#     if is_scalar:
#         return scale[0]
#     return scale


def read_sensitivity_fits(fits_file):
    table = QTable.read(fits_file)
    bin_edges = list(sorted(set(table['left_edge']) | set(table['right_edge']))) * u.TeV

    # log bin center
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_widths = np.diff(bin_edges)

    return table, bin_edges, bin_centers, bin_widths



def calculate_sensitivity(
        gammas,
        protons,
        min_energy,
        max_energy,
        gamma_prediction_cut=0.5,
        signal_region=0.01,
        target_spectrum=CrabSpectrum()
):

    if 'theta' not in gammas.columns:
        gammas['theta'] = coordinates.calculate_distance_theta(gammas, source_az=0 * u.deg, source_alt = 70 * u.deg).to(u.deg).value

    if 'theta' not in protons.columns:
        protons['theta'] = coordinates.calculate_distance_theta(protons, source_az=0 * u.deg, source_alt = 70 * u.deg).to(u.deg).value

    selected_gammas = gammas.query(f'gamma_prediction_mean >={gamma_prediction_cut}').copy()
    selected_protons = protons.query(f'gamma_prediction_mean >={gamma_prediction_cut}').copy()
    # print(selected_gammas.shape)
    # print(selected_protons.shape)

    # if len(selected_protons) < 50 or len(selected_gammas) < 50:
    #     return np.inf / (u.erg * u.s * u.cm**2)

    # on, off, alpha = coordinates.split_on_off(
    #     selected_protons,
    #     selected_gammas,
    #     on_region_radius=signal_region * u.deg
    # )
    #
    # n_off = off.weight.sum()
    # n_on = on.weight.sum()

    n_on, n_off = get_on_and_off_counts(
        selected_gammas,
        selected_protons,
        on_region_radius=signal_region
    )
    alpha = 1
    # if (n_off * alpha + 10 > n_on):
    #     return np.inf / (u.erg * u.s * u.cm**2)
    print(f'N_on = {n_on} N_off = {n_off}')
    relative_flux = relative_sensitivity(
        n_on,
        n_off,
        alpha=alpha,
    )
    # print(f'conf cut = {gamma_prediction_mean}, region = {signal_region},  relative_flux= {relative_flux}, lima = {li_ma_significance(n_on, n_off, alpha=alpha)}, aplphja= {alpha}')
    # print(f'n_on = {n_on}, n_off*alpha = {alpha*n_off}, n_off = {n_off}, ')

    # print(relative_flux)
    bin_center = np.sqrt(min_energy * max_energy) * u.TeV
    sens = target_spectrum.flux(bin_center) * relative_flux
    return sens


def get_on_and_off_counts(selected_gammas, selected_protons, on_region_radius):
    """ Get on and off counts from the signal region using a simpel theta**2 cut"""

    # estimate n_off by assuming that the background rate is constant within a
    # smallish theta area around 0. take the mean of the theta square histogram
    # to get a more stable estimate for n_off
    #
    # theta_square_cut = on_region_radius.to(u.deg).value**2
    #
    # H, _ = np.histogram(
    #     selected_protons.theta**2,
    #     bins=np.arange(0, 0.6, theta_square_cut),
    #     weights=selected_protons.weight
    # )
    # n_off = H.mean()
    n_off = selected_protons.query(f'theta < {on_region_radius.to(u.deg).value}')['weight'].sum()
    n_on = n_off + selected_gammas.query(f'theta < {on_region_radius.to(u.deg).value}')['weight'].sum()
    return n_on, n_off


def calculate_differential_sensitivity(
            gammas,
            protons,
            bin_edges,
            gamma_prediction_cut = 0.7,
            signal_region_radius = 0.25 * u.deg,
            target_spectrum=CrabSpectrum(),
):
    if 'theta' not in gammas.columns:
        gammas['theta'] = coordinates.calculate_distance_theta(gammas, source_az=0 * u.deg, source_alt=70 * u.deg).to(u.deg).value

    if 'theta' not in protons.columns:
        protons['theta'] = coordinates.calculate_distance_theta(protons, source_az=0 * u.deg, source_alt=70 * u.deg).to(u.deg).value

    gammas['energy_bin'] = pd.cut(gammas.mc_energy, bin_edges)
    protons['energy_bin'] = pd.cut(protons.mc_energy, bin_edges)

    rows = []
    for (bin, g,), ( _, p) in zip(gammas.groupby('energy_bin'), protons.groupby('energy_bin')):
        flux  = calculate_sensitivity(g, p, bin.left, bin.right, gamma_prediction_cut=gamma_prediction_cut, signal_region=signal_region_radius)

        d  = {'left_edge': bin.left, 'right_edge': bin.right, 'flux': flux.to(1 / (u.m**2 * u.s * u.TeV)).value}
        rows.append(d)

    t = Table(rows)
    t['flux'] = t['flux'] / (u.m**2 * u.s * u.TeV)
    t['radius'] = signal_region_radius
    t['threshold'] = gamma_prediction_cut
    t['left_edge'] = t['left_edge'] * u.TeV
    t['right_edge'] = t['right_edge'] * u.TeV
    return t

def optimize_differential_sensitivity(
            protons,
            gammas,
            bin_edges,
            target_spectrum=CrabSpectrum(),
            num_threads=-1,
):

    if 'theta' not in gammas.columns:
        gammas['theta'] = coordinates.calculate_distance_theta(gammas, source_az=0 * u.deg, source_alt=70 * u.deg).to(u.deg).value

    if 'theta' not in protons.columns:
        protons['theta'] = coordinates.calculate_distance_theta(protons, source_az=0 * u.deg, source_alt=70 * u.deg).to(u.deg).value

    gammas['energy_bin'] = pd.cut(gammas.mc_energy, bin_edges)
    protons['energy_bin'] = pd.cut(protons.mc_energy, bin_edges)

    if num_threads == -1:
        num_threads = multiprocessing.cpu_count() // 2

    args = [(gammas[gammas.energy_bin == bin], protons[protons.energy_bin == bin], bin) for bin in gammas.energy_bin.cat.categories]

    if num_threads > 1:
        with multiprocessing.Pool(processes=num_threads) as pool:
            results = pool.starmap(_find_best_sensitivity_in_bin, args)
    else:
        results = list(starmap(_find_best_sensitivity_in_bin, args))

    # multiply the whole thing by the proper unit. There must be a nicer way to do this.
    sensitivity = np.array([s[0].value for s in results]) * results[0][0].unit
    return sensitivity, _optimizer_result_to_table(results)


def _optimizer_result_to_table(result):
    d = [{'left_edge': bin.left, 'right_edge': bin.right, 'flux': flux.value, 'threshold': cut[0], 'radius': cut[1]} for flux, cut, bin in result]
    t = Table(d)
    t['flux'] = t['flux'] * result[0][0].unit
    t['radius'] = t['radius'] * u.deg
    t['left_edge'] = t['left_edge'] * u.TeV
    t['right_edge'] = t['right_edge'] * u.TeV
    return t


def _find_best_sensitivity_in_bin(g, p, bin):
    print(f'bin{bin}:{len(g)} gammas and {len(p)} protons ')
    min_energy, max_energy = bin.left, bin.right

    def f(x):
        return calculate_sensitivity(g, p, min_energy, max_energy, gamma_prediction_mean=x[0], signal_region=x[1]).value

    # ranges = (slice(0.0, 1, 0.025), slice(0.001, 0.08, 0.001))
    ranges = (slice(0.0, 1, 0.05))
    # Note: while it seems obviuous to use finish=optimize.fmin here. apparently it
    # tests invalid values. and then everything breaks. Negative theta cuts for
    # example
    res = optimize.brute(f, ranges, finish=None, full_output=True)

    cuts = res[0]
    print(f'best result in bin {bin} is: {cuts}')
    return calculate_sensitivity(g, p, min_energy, max_energy, gamma_prediction_mean=cuts[0], signal_region=cuts[1]), cuts, bin
