from scipy.optimize import minimize_scalar
from fact.analysis import li_ma_significance
import astropy.units as u
import numpy as np
import pandas as pd
from spectrum import CrabSpectrum
import multiprocessing
from joblib import delayed, Parallel
from scipy import optimize
import tqdm


@u.quantity_input(t_obs=u.hour, t_ref=u.hour)
def relative_sensitivity(
        n_on,
        n_off,
        alpha,
        target_significance=5,
):
    '''
    Calculate the relative sensitivity defined as the flux
    relative to the reference source that is detectable with
    significance in t_ref.

    Parameters
    ----------
    n_on: int or array-like
        Number of signal-like events for the on observations
    n_off: int or array-like
        Number of signal-like events for the off observations
    alpha: float
        Scaling factor between on and off observations.
        1 / number of off regions for wobble observations.
    target_significance: float
        Significance necessary for a detection

    Returns
    ----------
    The relative flux neccessary to detect the source with the given target significance.
    '''
    is_scalar = np.isscalar(n_on) and np.isscalar(n_off)

    if is_scalar:
        n_on = [n_on]
        n_off = [n_off]

    scale = []
    for on, off in zip(n_on, n_off):
        if on < off * alpha or off <= 10:
            scale.append(np.inf)
            continue

        def f(relative_flux):
            s = li_ma_significance((on - off) * relative_flux + off, off, alpha=alpha)
            return (target_significance - s)**2

        s = minimize_scalar(f, bounds=(1e-13, 300), method='bounded')

        scale.append(s.x)

    if is_scalar:
        return scale[0]
    return scale




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
        return np.inf / (u.erg * u.s * u.cm**2)

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
    """ Get on and off counts from the signal region using a simpel theta**2 cut"""

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


def find_differential_sensitivity(
            protons,
            gammas,
            bin_edges,
            target_spectrum=CrabSpectrum(),
            num_threads=-1
        ):

    gammas['energy_bin'] = pd.cut(gammas.mc_energy, bin_edges)
    protons['energy_bin'] = pd.cut(protons.mc_energy, bin_edges)

    if num_threads == -1:
        num_threads = multiprocessing.cpu_count() // 2

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
    return sensitivity


def find_best_sensitivity_in_bin(g, p):

    def f(x):
        return calculate_sensitivity(g, p, gamma_prediction_mean=x[0], signal_region=x[1]).value

    # ranges = (slice(0.0, 1, 0.025), slice(0.001, 0.08, 0.001))
    ranges = (slice(0.0, 1, 0.05), slice(0.001, 0.1, 0.0025))
    # Note: while it seems obviuous to use finish=optimize.fmin here. apparently it
    # tests invalid values. and then everything breaks. Negative theta cuts for
    # example
    res = optimize.brute(f, ranges, finish=None, full_output=True)

    cuts = res[0]
    return calculate_sensitivity(g, p, gamma_prediction_mean=cuts[0], signal_region=cuts[1])
