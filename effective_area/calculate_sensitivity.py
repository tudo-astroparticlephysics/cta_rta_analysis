import click
import astropy.units as u
from spectrum import CosmicRaySpectrum, CrabSpectrum, MCSpectrum, make_energy_bins
import fact.io
from sensitivity import optimize_differential_sensitivity, calculate_differential_sensitivity


@click.command()
@click.argument('gamma_input', type=click.Path(exists=True))
@click.argument('proton_input', type=click.Path(exists=True))
@click.argument('output', type=click.Path(exists=False))
@click.option('-n', '--n_bins', type=click.INT, default=20, help='energy bins to calculate')
@click.option('-j', '--n_jobs', type=click.INT, default=-1, help='number of threads to use in parallel')
@click.option('--optimize/--no-optimize', default=False, help='whether to optimize cuts per energy bin')
@click.option('-i', '--iterations', type=click.INT, default=1, help='number of iterations to perform for error calculation')
def main(
    gamma_input, proton_input, output,
    n_bins,
    n_jobs,
    optimize,
    iterations,
):
    '''
    Calculates a sensitivity curve vs real energy. For each energy bin it performs a gridsearch
    to find the theta and gamma_prediction_mean cuts that produce the highest sensitivity.
    '''

    t_obs = 50 * u.h
    e_min, e_max = 0.003 * u.TeV, 330 * u.TeV
    bin_edges, _, _ = make_energy_bins(e_min=e_min, e_max=e_max, bins=n_bins)

    columns = ['gamma_prediction_mean', 'az_prediction', 'alt_prediction', 'mc_alt', 'mc_az', 'mc_energy']

    gammas = fact.io.read_data(gamma_input, key='array_events', columns=columns)
    gammas = gammas.dropna()


    gamma_runs = fact.io.read_data(gamma_input, key='runs')
    mc_production_gamma = MCSpectrum.from_cta_runs(gamma_runs)

    protons = fact.io.read_data(proton_input, key='array_events', columns=columns)
    protons = protons.dropna()

    proton_runs = fact.io.read_data(proton_input, key='runs')
    mc_production_proton = MCSpectrum.from_cta_runs(proton_runs)

    crab = CrabSpectrum()
    cosmic = CosmicRaySpectrum()

    gammas['weight'] = mc_production_gamma.reweigh_to_other_spectrum(crab, gammas.mc_energy.values * u.TeV, t_assumed_obs=t_obs)
    protons['weight'] = mc_production_proton.reweigh_to_other_spectrum(cosmic, protons.mc_energy.values * u.TeV, t_assumed_obs=t_obs)


    if optimize:
        result_table = optimize_differential_sensitivity(gammas, protons, bin_edges=bin_edges, num_threads=n_jobs)
    else:
        result_table = calculate_differential_sensitivity(gammas, protons, bin_edges=bin_edges)

    result_table.write(output, overwrite=True)


if __name__ == '__main__':
    main()
