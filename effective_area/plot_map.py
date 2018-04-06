import click
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import fact.io
from astropy.coordinates import Angle
from coordinates import horizontal_to_skycoord, wrap_angles
from spectrum import MCSpectrum, CrabSpectrum, CosmicRaySpectrum


#

@click.command()
@click.argument('gammas_dl3', type=click.Path(exists=True))
@click.argument('protons_dl3', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
def main(gammas_dl3, protons_dl3, output):

    # t_obs = 3.6 * u.h
    cut = 0.0

    gammas = fact.io.read_data(gammas_dl3, key='array_events')
    gammas = gammas.query(f'gamma_prediction_mean > {cut}')
    #
    # gamma_runs = fact.io.read_data(gammas_dl3, key='runs')
    # mc_production_gamma = MCSpectrum.from_cta_runs(gamma_runs)

    protons = fact.io.read_data(protons_dl3, key='array_events')
    protons = protons.query(f'gamma_prediction_mean > {cut}')

    # print(f'Plotting {len(protons)} protons and {len(gammas)} gammas.')
    # proton_runs = fact.io.read_data(protons_dl3, key='runs')
    # mc_production_proton = MCSpectrum.from_cta_runs(proton_runs)
    #
    # crab = CrabSpectrum()
    # cosmic = CosmicRaySpectrum()

    # gamma_weights = mc_production_gamma.reweigh_to_other_spectrum(crab, gammas.mc_energy.values * u.TeV, t_assumed_obs=t_obs)
    # proton_weights = mc_production_proton.reweigh_to_other_spectrum(cosmic, protons.mc_energy.values * u.TeV, t_assumed_obs=t_obs)

    # az = gammas.az_prediction.values * u.deg
    # alt = gammas.alt_prediction.values * u.deg
    # c_gamma = horizontal_to_skycoord(alt, az)

    fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(10, 5))

    az = protons.mc_az.values * u.rad
    alt = protons.mc_alt.values * u.rad
    c = horizontal_to_skycoord(alt, az)

    ra = c.icrs.ra
    dec = c.icrs.dec
    # ax1.scatter(ra.degree, dec.degree, label='proton', s=2, alpha=0.5)
    ax1.scatter(*wrap_angles(alt, az), label=f'proton ({len(az)})', s=2, alpha=0.5)

    az = gammas.mc_az.values * u.rad
    alt = gammas.mc_alt.values * u.rad
    c = horizontal_to_skycoord(alt, az)

    ra = c.icrs.ra
    dec = c.icrs.dec
    # ax1.scatter(ra.degree, dec.degree, label='gamma', s=2, alpha=0.5)
    ax1.scatter(*wrap_angles(alt, az), label=f'gamma ({len(az)})', s=8, alpha=0.5)
    ax1.set_title('True Direction')


    az = protons.az_prediction.values * u.rad
    alt = protons.alt_prediction.values * u.rad
    c = horizontal_to_skycoord(alt, az)

    ra = c.icrs.ra
    dec = c.icrs.dec
    # ax2.scatter(ra.degree, dec.degree, label='proton', s=2, alpha=0.5)
    ax2.scatter(*wrap_angles(alt, az), label=f'proton ({len(az)})', s=2, alpha=0.5)

    az = gammas.az_prediction.values * u.rad
    alt = gammas.alt_prediction.values * u.rad

    c = horizontal_to_skycoord(alt, az)

    ra = c.icrs.ra
    dec = c.icrs.dec
    # ax2.scatter(ra.degree, dec.degree, label='gamma', s=2, alpha=0.5)
    ax2.scatter(*wrap_angles(alt, az), label=f'gamma ({len(az)})', s=2, alpha=0.5)
    ax2.set_title('Reconstruncted Direction')
    #
    ax1.set_ylim([-10, 10])
    ax1.set_xlim([60, 80])
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_ylim(ax1.get_ylim())
    # az = protons.az_prediction.values * u.deg
    # alt = protons.alt_prediction.values * u.deg
    # c_proton = horizontal_to_skycoord(alt, az)
    #
    # ra = c_proton.icrs.ra
    # dec = c_proton.icrs.dec
    # plt.scatter(ra.degree, dec.degree, label='proton_prediction', s=1, alpha=0.5)
    ax2.set_xlabel('alt')
    ax2.set_ylabel('az')

    ax1.set_xlabel('alt')
    ax1.set_ylabel('az')


    ax1.legend()
    ax2.legend()

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
