import click
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import fact.io
import coordinates
from spectrum import MCSpectrum, CrabSpectrum, CosmicRaySpectrum


@click.command()
@click.argument('gammas_dl3', type=click.Path(exists=True))
@click.argument('protons_dl3', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
def main(gammas_dl3, protons_dl3, output):

    t_obs = 0.5 * u.h
    cut = 0.0

    gammas = fact.io.read_data(gammas_dl3, key='array_events')
    gammas = gammas.dropna()


    gamma_runs = fact.io.read_data(gammas_dl3, key='runs')
    mc_production_gamma = MCSpectrum.from_cta_runs(gamma_runs)

    protons = fact.io.read_data(protons_dl3, key='array_events')
    protons = protons.dropna()

    # print(f'Plotting {len(protons)} protons and {len(gammas)} gammas.')
    proton_runs = fact.io.read_data(protons_dl3, key='runs')
    mc_production_proton = MCSpectrum.from_cta_runs(proton_runs)


    crab = CrabSpectrum()
    cosmic = CosmicRaySpectrum()


    gammas['weight'] = mc_production_gamma.reweigh_to_other_spectrum(crab, gammas.mc_energy.values * u.TeV, t_assumed_obs=t_obs)
    protons['weight'] = mc_production_proton.reweigh_to_other_spectrum(cosmic, protons.mc_energy.values * u.TeV, t_assumed_obs=t_obs)

    # min = gammas.mc_energy.min() * u.TeV
    # max = gammas.mc_energy.max() * u.TeV
    # n_crab = crab.expected_events(min, max, area=mc_production_gamma.generation_area, t_obs=t_obs)

    # min = protons.mc_energy.min() * u.TeV
    # max = protons.mc_energy.max() * u.TeV
    # n_cosmic = cosmic.expected_events(min, max, area=mc_production_proton.generation_area, t_obs=t_obs, solid_angle=mc_production_proton.generator_solid_angle)


    gammas_gammalike = gammas.query(f'gamma_prediction_mean > {cut}')
    protons_gammalike = protons.query(f'gamma_prediction_mean > {cut}')


    gammas_gammalike['theta'] = coordinates.calculate_distance_theta(gammas_gammalike)
    protons_gammalike['theta'] = coordinates.calculate_distance_theta(protons_gammalike)

    bins = np.linspace(0, 0.3, 20)
    kwargs = {'histtype': 'step', 'lw': 2.0}


    plt.hist(gammas_gammalike['theta']**2, bins=bins, label='gamma', weights=gammas_gammalike.weight, **kwargs)
    h, _, _ = plt.hist(protons_gammalike['theta']**2, bins=bins, label='proton', weights=protons_gammalike.weight, **kwargs)
    plt.suptitle(f'Observation Time {t_obs} \n Prediction threshold {cut}')
    print(f'Mean Bkg per bin {h.mean()}, {h.std()}')

    if output:
        plt.savefig(output)
    else:
        plt.show()

#
#
#
# @click.command()
# @click.argument('predicted_gammas', type=click.Path(exists=True, dir_okay=False,))
# @click.argument('predicted_protons', type=click.Path(exists=True, dir_okay=False,))
# @click.argument('outputfile', type=click.Path(exists=False, dir_okay=False,))
# @click.option('-n', '--n_bins', type=click.INT, default=30, help='theta bin')
# def main(
#     predicted_gammas,
#     predicted_protons,
#     mc_production_information,
#     outputfile,
#     n_bins,
#     sample_fraction,
# ):
#     '''
#     Plot the famous theta square curve.
#     '''
#     t_obs = 3.6 * u.h
#
#     # read the gammas and weight them accroding to the crab spectrum
#     gammas, N, e_min, e_max, area = cta_io.read_events(predicted_gammas, mc_production_information)
#
#     mc_gamma = power_law.MCSpectrum(
#         e_min=e_min,
#         e_max=e_max,
#         total_showers_simulated=N * sample_fraction,
#         generation_area=area,
#     )
#
#     crab = power_law.CrabSpectrum()
#     energies = gammas.energy.values * u.TeV
#
#     gammas['weight'] = crab.weight(
#         energies,
#         mc_spectrum=mc_gamma,
#         t_assumed_obs=t_obs,
#     )
#
#     # read the protons and weight them accroding to the cosmic ray spectrum
#     protons, N, e_min, e_max, area = cta_io.read_events(predicted_protons, mc_production_information)
#     cosmic_spectrum = power_law.CosmicRaySpectrum()
#     mc_proton = power_law.MCSpectrum(
#         e_min=e_min,
#         e_max=e_max,
#         total_showers_simulated=N,
#         generation_area=area,
#         generator_solid_angle=6 * u.deg
#     )
#
#     energies = protons.energy.values*u.TeV
#
#     protons['weight'] = cosmic_spectrum.weight(
#         energies,
#         mc_spectrum=mc_proton,
#         t_assumed_obs=t_obs,
#     )
#
#     # select gamma-like events from both samples and plot the theta histogram
#     selected_protons = protons.query('gammaness >= 0.7')
#     selected_gammas = gammas.query('gammaness >= 0.7')
#
#     _, edges, _ = plt.hist(
#                     selected_protons.theta_deg**2,
#                     bins=n_bins,
#                     range=[0, 0.2],
#                     weights=selected_protons.weight,
#                     histtype='step',
#                 )
#     plt.hist(
#         selected_gammas.theta_deg**2,
#         bins=edges,
#         weights=selected_gammas.weight,
#         histtype='step',
#     )
#     plt.xlabel('$(\mathrm{Theta} / \mathrm{degree})^2$')
#     plt.ylabel('Expected events in {} '.format(t_obs))
#     plt.savefig(outputfile)
#
#
if __name__ == '__main__':
    main()
