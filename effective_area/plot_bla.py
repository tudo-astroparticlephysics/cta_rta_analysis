import click
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import coordinates
import fact.io


@click.command()
@click.argument('protons_dl3', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
def main(protons_dl3, output):
    #
    df = fact.io.read_data(protons_dl3, key='array_events')
    df = df.dropna()
    #
    # # print(f'Plotting {len(protons)} protons and {len(gammas)} gammas.')
    # proton_runs = fact.io.read_data(protons_dl3, key='runs')
    # alt_pointing, az_pointing = 70, 0
    #
    # az = Angle(df.az_prediction.values, unit=u.rad).wrap_at(180*u.deg).degree
    # alt = Angle(df.alt_prediction.values, unit=u.rad).degree



    # protons['weight'] = mc_production_proton.reweigh_to_other_spectrum(cosmic, protons.mc_energy.values * u.TeV, t_assumed_obs=t_obs)

    # protons_gammalike = protons.query(f'gamma_prediction_mean > {cut}')
    # protons_gammalike = protons


    distance = coordinates.calculate_distance_theta(df, source_az=0 * u.deg, source_alt=70 * u.deg).to('deg').value

    # distance = np.sqrt((az_pointing - az)**2 + (alt_pointing - alt)**2)

    bins = np.arange(0, 13, 0.5)
    fig, [ax1, ax2] = plt.subplots(2, 1)
    ax1.hist(distance, bins=bins, label='proton',)

    bins = np.arange(0, 13, 0.5)
    ax2.hist(distance**2, bins=bins, label='proton',)

    #
    # bins = np.arange(0, (protons_gammalike['theta']**2).max(), 0.01)
    # h, _ = np.histogram(protons_gammalike['theta']**2, bins=bins,)
    #
    # print(f'Mean Bkg per bin {h.mean()}, {h.std()}')
    # # plt.axhline(h.mean(), color='gray')

    # plt.suptitle(f'Observation Time {t_ob/s} \n Prediction threshold {cut}')
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
