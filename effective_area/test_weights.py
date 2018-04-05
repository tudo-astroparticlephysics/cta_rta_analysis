import click
import numpy as np
import astropy.units as u
import fact.io
from spectrum import MCSpectrum, CrabSpectrum, CosmicRaySpectrum



@click.command()
@click.argument('input_gammas', type=click.Path(exists=True))
@click.argument('input_protons', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
def main(input_gammas, input_protons, output):

    t_obs = 3.6 * u.h
    #
    # runs = fact.io.read_data(input_gammas, key='runs')
    # mc_production_gamma = MCSpectrum.from_cta_runs(runs)
    #
    # crab = CrabSpectrum()
    #
    # energies = mc_production_gamma.draw_energy_distribution(0.1 * u.TeV, 100 * u.TeV, size=10000)
    # w1 = crab.weight(energies, mc_production_gamma, t_assumed_obs=t_obs)

    # now test the protons
    runs = fact.io.read_data(input_protons, key='runs')
    mc_production_proton = MCSpectrum.from_cta_runs(runs)
    energies = mc_production_proton.draw_energy_distribution(0.1 * u.TeV, 100 * u.TeV, size=10000)
    
    cosmic = CosmicRaySpectrum()

    w1 = cosmic.weight(energies, mc_production_proton, t_assumed_obs=t_obs)



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
