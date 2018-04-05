import click
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import fact.io
from astropy.coordinates import Angle
from spectrum import MCSpectrum, CrabSpectrum, CosmicRaySpectrum



def calculate_distance(df):
    alt = Angle(df.alt_prediction.values * u.rad).degree
    mc_alt = Angle(df.mc_alt.values * u.rad).degree

    az = Angle(df.az_prediction.values * u.rad).wrap_at(180 * u.deg).degree
    mc_az = Angle(df.mc_az.values * u.rad).wrap_at(180 * u.deg).degree

    return np.sqrt((alt - mc_alt)**2 + (az - mc_az)**2)


@click.command()
@click.argument('gammas_dl3', type=click.Path(exists=True))
@click.argument('protons_dl3', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
def main(gammas_dl3, protons_dl3, output):

    t_obs = 3.6 * u.h
    cut = 0.0

    gammas = fact.io.read_data(gammas_dl3, key='array_events')
    gammas = gammas.query(f'gamma_prediction_mean > {cut}')

    gamma_runs = fact.io.read_data(gammas_dl3, key='runs')
    mc_production_gamma = MCSpectrum.from_cta_runs(gamma_runs)

    protons = fact.io.read_data(protons_dl3, key='array_events')
    protons = protons.query(f'gamma_prediction_mean > {cut}')

    print(f'Plotting {len(protons)} protons and {len(gammas)} gammas.')
    proton_runs = fact.io.read_data(protons_dl3, key='runs')
    mc_production_proton = MCSpectrum.from_cta_runs(proton_runs)

    crab = CrabSpectrum()
    cosmic = CosmicRaySpectrum()

    gamma_weights = mc_production_gamma.reweigh_to_other_spectrum(crab, gammas.mc_energy.values * u.TeV, t_assumed_obs=t_obs)
    proton_weights = mc_production_proton.reweigh_to_other_spectrum(cosmic, protons.mc_energy.values * u.TeV, t_assumed_obs=t_obs)

    gammas['theta'] = calculate_distance(gammas)
    protons['theta'] = calculate_distance(protons)

    bins = np.linspace(0, 0.2, 40)
    kwargs = {'histtype': 'step', 'lw': 2.0}


    plt.hist(gammas['theta']**2, bins=bins, label='gamma', **kwargs)
    plt.hist(protons['theta']**2, bins=bins, label='proton', **kwargs)
    plt.legend()

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
