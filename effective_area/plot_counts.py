import click
import matplotlib.pyplot as plt
import fact.io
import astropy.units as u
import coordinates
from spectrum import make_energy_bins
from astropy.coordinates import Angle, SkyCoord
from regions import CircleSkyRegion
import numpy as np
from spectrum import CosmicRaySpectrum, CrabSpectrum, MCSpectrum
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages


def calculate_distance(df):
    alt = Angle(df.alt_prediction.values * u.rad).degree
    mc_alt = Angle(df.mc_alt.values * u.rad).degree

    az = Angle(df.az_prediction.values * u.rad).wrap_at(180 * u.deg).degree
    mc_az = Angle(df.mc_az.values * u.rad).wrap_at(180 * u.deg).degree

    return np.sqrt((alt - mc_alt)**2 + (az - mc_az)**2)



@click.command()
@click.argument('gammas_dl3', type=click.Path(exists=True))
@click.argument('protons_dl3', type=click.Path(exists=True))
@click.argument('output_pdf', type=click.Path(exists=False))
@click.option('-b', '--bins', default=10, show_default=True)
@click.option('-t', '--threshold', default=0.5, show_default=True, help='prediction threshold')
def main(gammas_dl3, protons_dl3, output_pdf, bins, threshold):
    cols = ['az_prediction', 'alt_prediction', 'mc_energy']


    gammas = fact.io.read_data(gammas_dl3, key='array_events')
    print(f'Reading {len(gammas)} gammas')
    # import IPython; IPython.embed()

    gammas = gammas.dropna(subset=cols)

    protons = fact.io.read_data(protons_dl3, key='array_events')
    print(f'Reading {len(protons)} protons')
    protons = protons.dropna(subset=cols)

    gamma_runs = fact.io.read_data(gammas_dl3, key='runs')
    mc_production_gamma = MCSpectrum.from_cta_runs(gamma_runs)

    proton_runs = fact.io.read_data(protons_dl3, key='runs')
    mc_production_proton = MCSpectrum.from_cta_runs(proton_runs)

    crab = CrabSpectrum()
    cosmic = CosmicRaySpectrum()
    t_obs = 60 * u.s
    gammas['weight'] = mc_production_gamma.reweigh_to_other_spectrum(crab, gammas.mc_energy.values * u.TeV, t_assumed_obs=t_obs)
    protons['weight'] = mc_production_proton.reweigh_to_other_spectrum(cosmic, protons.mc_energy.values * u.TeV, t_assumed_obs=t_obs)

    with PdfPages(output_pdf) as pdf:

        plt.figure()
        bin_edges, _, _ = make_energy_bins(gammas.mc_energy.values * u.TeV, bins=20)
        plt.hist(gammas.mc_energy, bins=bin_edges, histtype='step', lw=2, label='gammas true energy')
        plt.hist(protons.mc_energy, bins=bin_edges, histtype='step', lw=2, label='proton true energy')
        plt.legend()

        plt.xscale('log')
        plt.xlabel(r'$Energy /  \mathrm{TeV}$')
        plt.ylabel('Counts')
        pdf.savefig()
        plt.close()

        plt.figure()
        gammas['theta'] = calculate_distance(gammas)
        protons['theta'] = calculate_distance(protons)

        bins = np.linspace(0, 25.2, 40)
        kwargs = {'histtype': 'step', 'lw': 2.0}
        plt.hist(gammas['theta'], bins=bins, label='gamma', density=True, **kwargs)
        plt.hist(protons['theta'], bins=bins, label='proton', density=True, **kwargs)
        plt.legend()
        plt.xlabel(r'$Distance /  Degree$')
        plt.ylabel('Normalized Counts')
        pdf.savefig()
        plt.close()


        plt.figure()
        plt.hist(gammas['theta']**2, bins=bins, label='gamma', density=True, **kwargs)
        plt.hist(protons['theta']**2, bins=bins, label='proton', density=True, **kwargs)
        plt.legend()
        plt.xlabel(r'$Distance Squared /  Degree^2$')
        plt.ylabel('Normalized Counts')
        pdf.savefig()
        plt.close()

        plt.figure()
        bins = np.linspace(0, 0.30, 40)
        plt.hist(gammas['theta']**2, bins=bins, label='gamma', density=True, **kwargs)
        plt.hist(protons['theta']**2, bins=bins, label='proton', density=True, **kwargs)
        plt.legend()
        plt.xlabel(r'$Distance Squared /  Degree^2$')
        plt.ylabel('Normalized Counts')
        pdf.savefig()
        plt.close()

        plt.figure()
        bins = np.linspace(0, 0.30, 40)
        on = gammas.append(protons)
        plt.hist(gammas['theta']**2, bins=bins, label='gammas', weights=gammas.weight, **kwargs)
        plt.hist(protons['theta']**2, bins=bins, label='proton', weights=protons.weight, **kwargs)
        plt.legend()
        plt.xlabel(r'$Distance Squared /  Degree^2$')
        plt.ylabel('weighted Counts')
        pdf.savefig()
        plt.close()

        plt.figure()
        bins = np.linspace(0, 0.30, 40)
        on = gammas.append(protons)
        plt.hist(on['theta']**2, bins=bins, label='on (gammas + protons)', weights=on.weight, **kwargs)
        plt.hist(protons['theta']**2, bins=bins, label='proton', weights=protons.weight, **kwargs)
        plt.legend()
        plt.xlabel(r'$Distance Squared /  Degree^2$')
        plt.ylabel('weighted Counts')
        pdf.savefig()
        plt.close()

        fig, ax = plt.subplots()

        c = coordinates.skyccords_from_dl3_table(gammas).icrs
        ra = c.ra.deg
        dec = c.dec.deg

        ra_bins = np.linspace(11, 17, 50)
        dec_bins = np.linspace(41, 47, 50)

        _, _, _, im = ax.hist2d(ra, dec, bins=[ra_bins, dec_bins], cmap='magma')
        ax.set_aspect('equal',)
        plt.colorbar(im, ax=ax)
        plt.title('gammas')
        plt.xlabel(r'$Right Ascension /  Degree$')
        plt.ylabel(r'$Declination /  Degree$')
        pdf.savefig()
        plt.close()



        fig, ax = plt.subplots()

        c = coordinates.skyccords_from_dl3_table(protons).icrs
        ra = c.ra.deg
        dec = c.dec.deg

        _, _, _, im = ax.hist2d(ra, dec, bins=[ra_bins, dec_bins], cmap='magma')
        ax.set_aspect('equal')
        plt.title('protons')
        plt.colorbar(im, ax=ax)
        plt.xlabel(r'$Right Ascension /  Degree$')
        plt.ylabel(r'$Declination /  Degree$')
        pdf.savefig()
        plt.close()



        fig, ax = plt.subplots()
        on, off, _ = coordinates.split_on_off(gammas, protons)
        c = coordinates.skyccords_from_dl3_table(on).icrs
        _, _, _, im = ax.hist2d(c.ra.deg, c.dec.deg, bins=[ra_bins, dec_bins], cmap='magma')
        ax.set_aspect('equal')
        plt.title('on region')
        plt.colorbar(im, ax=ax)
        pdf.savefig()
        plt.close()


        fig, ax = plt.subplots()
        c = coordinates.skyccords_from_dl3_table(off).icrs
        _, _, _, im = ax.hist2d(c.ra.deg, c.dec.deg, bins=[ra_bins, dec_bins], cmap='magma')
        ax.set_aspect('equal')
        plt.colorbar(im, ax=ax)
        plt.title('off region')
        pdf.savefig()
        plt.close()


        plt.figure()
        plt.hist(on.mc_energy, bins=bin_edges, histtype='step', lw=2, label='on region true energy', density=True)
        plt.hist(off.mc_energy, bins=bin_edges, histtype='step', lw=2, label='off region true energy',density=True)

        plt.hist(gammas.mc_energy, bins=bin_edges, histtype='step', lw=2, label='gammas true energy', color='gray', density=True)
        plt.hist(protons.mc_energy, bins=bin_edges, histtype='step', lw=2, label='proton true energy', color='lightgray', density=True)
        plt.legend()

        plt.xscale('log')
        plt.xlabel(r'$Energy /  \mathrm{TeV}$')
        plt.ylabel(' Normalized  Counts')
        pdf.savefig()
        plt.close()


        gammas = gammas.query(f'gamma_prediction_mean > {threshold}')
        protons = protons.query(f'gamma_prediction_mean > {threshold}')

        plt.figure()
        bins = np.linspace(0, 0.22, 40)
        plt.hist(gammas['theta']**2, bins=bins, label='gamma', density=True, **kwargs)
        plt.hist(protons['theta']**2, bins=bins, label='proton', density=True, **kwargs)
        plt.legend()
        plt.xlabel(r'$Distance Squared /  Degree^2$')
        plt.ylabel(f'Normalized Counts t > {threshold}')
        pdf.savefig()
        plt.close()

        fig, ax = plt.subplots()
        on, off, _ = coordinates.split_on_off(gammas, protons)
        c = coordinates.skyccords_from_dl3_table(on).icrs
        _, _, _, im = ax.hist2d(c.ra.deg, c.dec.deg, bins=[ra_bins, dec_bins], cmap='magma')
        ax.set_aspect('equal')
        plt.colorbar(im, ax=ax)
        plt.title(f'{len(on)} gammalike events in on region threshold={threshold}')
        pdf.savefig()
        plt.close()


        fig, ax = plt.subplots()
        c = coordinates.skyccords_from_dl3_table(off).icrs
        _, _, _, im = ax.hist2d(c.ra.deg, c.dec.deg, bins=[ra_bins, dec_bins], cmap='magma')
        ax.set_aspect('equal')
        plt.colorbar(im, ax=ax)
        plt.title(f'{len(off)} gammalike events in off region threshold={threshold}')
        pdf.savefig()
        plt.close()


        plt.figure()
        plt.title('gammalike events')
        plt.hist(on.mc_energy, bins=bin_edges, histtype='step', lw=2, label='on region true energy', density=True)
        plt.hist(off.mc_energy, bins=bin_edges, histtype='step', lw=2, label='off region true energy', density=True)

        plt.hist(gammas.mc_energy, bins=bin_edges, histtype='step', lw=2, label='gammas true energy', color='gray', density=True)
        plt.hist(protons.mc_energy, bins=bin_edges, histtype='step', lw=2, label='proton true energy', color='lightgray', density=True)
        plt.legend()

        plt.xscale('log')
        plt.xlabel(r'$Energy /  \mathrm{TeV}$')
        plt.ylabel('Normalized Counts')
        pdf.savefig()
        plt.close()



        plt.figure()
        plt.title('gammalike events')
        plt.hist(on.mc_energy, bins=bin_edges, histtype='step', lw=2, label='on region true energy', )
        plt.hist(off.mc_energy, bins=bin_edges, histtype='step', lw=2, label='off region true energy',)

        plt.hist(gammas.mc_energy, bins=bin_edges, histtype='step', lw=2, label='gammas true energy', color='gray', )
        plt.hist(protons.mc_energy, bins=bin_edges, histtype='step', lw=2, label='proton true energy', color='lightgray', )
        plt.legend()

        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$Energy /  \mathrm{TeV}$')
        plt.ylabel('Counts')
        pdf.savefig()
        plt.close()



if __name__ == "__main__":
    main()
