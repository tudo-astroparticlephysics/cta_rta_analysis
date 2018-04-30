import click
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import fact.io
from fact.analysis import li_ma_significance
import coordinates
from spectrum import MCSpectrum, CrabSpectrum, CosmicRaySpectrum


@click.command()
@click.argument('gammas_dl3', type=click.Path(exists=True))
@click.argument('protons_dl3', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
def main(gammas_dl3, protons_dl3, output):

    t_obs = 5 * u.h
    prediction_cut = 0.55
    signal_region = 0.185

    # cut = 0.85

    gammas = fact.io.read_data(gammas_dl3, key='array_events')
    gammas = gammas.dropna()


    gamma_runs = fact.io.read_data(gammas_dl3, key='runs')
    mc_production_gamma = MCSpectrum.from_cta_runs(gamma_runs)

    protons = fact.io.read_data(protons_dl3, key='array_events')
    print(len(protons))
    protons = protons.dropna()

    # print(f'Plotting {len(protons)} protons and {len(gammas)} gammas.')
    proton_runs = fact.io.read_data(protons_dl3, key='runs')
    print(len(proton_runs))
    mc_production_proton = MCSpectrum.from_cta_runs(proton_runs)


    crab = CrabSpectrum()
    cosmic = CosmicRaySpectrum()


    gammas['weight'] = mc_production_gamma.reweigh_to_other_spectrum(crab, gammas.mc_energy.values * u.TeV, t_assumed_obs=t_obs)
    protons['weight'] = mc_production_proton.reweigh_to_other_spectrum(cosmic, protons.mc_energy.values * u.TeV, t_assumed_obs=t_obs)

    # significance = 0
    # prediction_cut = 0
    # for cut in np.linspace(0, 1, 20):
    #
    #     gammas_gammalike = gammas.query(f'gamma_prediction_mean > {cut}').copy()
    #     protons_gammalike = protons.query(f'gamma_prediction_mean > {cut}').copy()
    #     # print(f'{len(gammas_gammalike)} gammas affter cut and {len(protons_gammalike)} after cut')
    #
    #
    #     gammas_gammalike['theta'] = coordinates.calculate_distance_theta(gammas_gammalike, source_az=0 * u.deg, source_alt=70 * u.deg).to(u.deg).value
    #     protons_gammalike['theta'] = coordinates.calculate_distance_theta(protons_gammalike, source_az=0 * u.deg, source_alt=70 * u.deg).to(u.deg).value
    #
    #     bin_width = 0.1**2
    #     bins = np.arange(0, 0.3, bin_width)
    #     kwargs = {'histtype': 'step', 'lw': 2.0}
    #
    #     on = gammas_gammalike.append(protons_gammalike)
    #     off = protons_gammalike
    #
    #     bins = np.arange(0, 0.3, signal_region**2)
    #     h, _ = np.histogram(off['theta']**2, bins=bins, weights=off.weight)
    #     # print(f'Mean Bkg {h.mean()}, {h.std()}')
    #
    #     n_off = h.mean()
    #     n_on = on.query(f'theta**2 <= {signal_region**2}').weight.sum()
    #     s = li_ma_significance(n_on, n_off, alpha=1)
    #     print(s, cut)
    #     if s > significance:
    #         prediction_cut = cut
    #

    gammas_gammalike = gammas.query(f'gamma_prediction_mean > {prediction_cut}').copy()
    protons_gammalike = protons.query(f'gamma_prediction_mean > {prediction_cut}').copy()
    # print(f'{len(gammas_gammalike)} gammas affter cut and {len(protons_gammalike)} after cut')


    gammas_gammalike['theta'] = coordinates.calculate_distance_theta(gammas_gammalike, source_az=0 * u.deg, source_alt=70 * u.deg).to(u.deg).value
    protons_gammalike['theta'] = coordinates.calculate_distance_theta(protons_gammalike, source_az=0 * u.deg, source_alt=70 * u.deg).to(u.deg).value

    bin_width = 0.1**2
    bins = np.arange(0, 0.3, bin_width)
    kwargs = {'histtype': 'step', 'lw': 2.0}

    on = gammas_gammalike.append(protons_gammalike)
    off = protons_gammalike
    plt.hist(on['theta']**2, bins=bins, label='on events', weights=on.weight, **kwargs)
    plt.hist(off['theta']**2, bins=bins, label='off events', weights=off.weight, **kwargs)
    plt.axvline(signal_region**2, color='gray')



    # print(f'N_on  = {n_on},        N_off  = {n_off} ')
    # print(f'Observation Time {t_obs} \n Prediction threshold {cut} \n Significance {significance}')

    # plt.suptitle(f'Observation Time {t_obs} \n Prediction threshold {cut} \n Significance {significance}')
    if output:
        plt.savefig(output)
    else:
        plt.show()


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
