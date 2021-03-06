import numpy as np
import click
import matplotlib.pyplot as plt
import astropy.units as u
from sensitivity import read_sensitivity_fits
from spectrum import CrabSpectrum
import pandas as pd



@u.quantity_input(bin_edges=u.TeV,)
def plot_sensitivity(bin_edges, sensitivity, sensitivity_std=None, ax=None, scale=True, **kwargs):

    y_error = None
    if not ax:
        ax = plt.gca()

    bin_center = np.sqrt(bin_edges[:-1] * bin_edges[1:])
    xerr = [np.abs(bin_edges[:-1] - bin_center).value, np.abs(bin_edges[1:] - bin_center).value]

    # sensitivity = sensitivity.to(1 / (u.erg * u.s * u.cm**2))
    # import IPython; IPython.embed()

    if sensitivity_std is not None:
        y_error = True
        # sensitivity_std = sensitivity_std.to(1 / (u.erg * u.s * u.cm**2))

    if scale:
        sensitivity = sensitivity * bin_center.to('erg')**2
        if y_error:
            sensitivity_std = sensitivity_std * bin_center.to('erg')**2

    ax.errorbar(
        bin_center.value,
        sensitivity.value,
        xerr=xerr,
        yerr=sensitivity_std.value / 2 if y_error else None,
        marker='.',
        linestyle='',
        capsize=0,
        **kwargs,
    )

    return ax


@u.quantity_input(e_min=u.TeV, e_max=u.TeV)
def plot_spectrum(spectrum, e_min, e_max, ax=None, scale=True, **kwargs):

    if not ax:
        _, ax = plt.subplots(1)

    e = np.linspace(e_min, e_max, 1000)
    flux = spectrum.flux(e).to(1 / (u.erg * u.s * u.cm**2))

    if scale:
        flux = flux * e.to('erg')**2

    ax.plot(
        e,
        flux,
        linestyle='--',
        **kwargs,
    )
    return ax




@click.command()
@click.argument('input_fits', type=click.Path(exists=True), nargs=-1)
@click.option('-l', '--label', multiple=True)
@click.option('-c', '--color', multiple=True)
@click.option('-o', '--output', type=click.Path(exists=False))
@click.option('--reference/--no-reference', default=True)
def main(
    input_fits, label, color,
    output, reference
):
    '''
    Plots a sensitivity curve vs real energy. For each energy bin it performs a gridsearch
    to find the theta and gamma_prediction_mean cuts that produce the highest sensitivity.
    '''


    if label and len(input_fits) != len(label):
        print('Must pass as many labels as gamma files as proton files')

    if color and label and len(label) != len(color):
        print('Must pass as many colors as labels')

    t_obs = 50 * u.h

    if not label:
        label = input_fits

    if not color:
        color = [None] * len(input_fits)

    for fits_file, l, c in zip(input_fits, label, color):
        table, bin_edges, _, _ = read_sensitivity_fits(fits_file)

        sens = table['flux'].to(1 / (u.erg * u.s * u.cm**2))
        if 'flux_std' in table.colnames:
            flux_std = table['flux_std'].to(1 / (u.erg * u.s * u.cm**2))
        else:
            flux_std = None

        ax = plot_sensitivity(bin_edges, sens, sensitivity_std=flux_std, ax=None, label=l, color=c)

    plot_spectrum(CrabSpectrum(), 0.003 * u.TeV, 330 * u.TeV, ax=ax, color='gray')

    if reference:
        path = 'resources/ascii/CTA-Performance-prod3b-v1-South-20deg-50h-DiffSens.txt'
        df = pd.read_csv(path, delimiter='\t\t', skiprows=10, names=['e_min', 'e_max', 'sensitivity'], engine='python')
        bin_edges = sorted(list(set(df.e_min) | set(df.e_max))) * u.TeV
        sensitivity = df.sensitivity.values * u.erg/(u.cm**2 * u.s)
        ax = plot_sensitivity(bin_edges, sensitivity, ax=ax, scale=False, label='prod3B reference', color='black')


    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel(r'$ \mathrm{photons} / \mathrm{erg s} \mathrm{cm}^2$ in ' + str(t_obs.to('h')) + ' hours' )
    ax.set_ylabel(r'$ E^2 \cdot \mathrm{photons} \quad \mathrm{erg} /( \mathrm{s} \quad  \mathrm{cm}^2$ )  in ' + str(t_obs.to('h')) )
    ax.set_xlabel(r'$E /  \mathrm{TeV}$')


    if label[0]:
        plt.legend()

    if output:
        plt.savefig(output)
    else:
        plt.show()



if __name__ == '__main__':
    main()
