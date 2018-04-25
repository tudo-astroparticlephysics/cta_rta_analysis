import click
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
from scipy.stats import binned_statistic
import fact.io


@click.command()
@click.argument('input_dl3_file', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
@click.option('-t', '--title', default=None)
@click.option('-c', '--color', default='#0a4f88')
def main(input_dl3_file, output, title, color):
    df = fact.io.read_data(input_dl3_file, key='array_events').dropna()

    alt = Angle(df.alt_prediction.values * u.rad).degree
    mc_alt = Angle(df.mc_alt.values * u.rad).degree

    az = Angle(df.az_prediction.values * u.rad).wrap_at(180 * u.deg).degree
    mc_az = Angle(df.mc_az.values * u.rad).wrap_at(180 * u.deg).degree

    distance = np.sqrt((alt - mc_alt)**2 + (az - mc_az)**2)

    bins = np.logspace(np.log10(df.mc_energy.min()), np.log10(df.mc_energy.max()), 10)
    df['energy_bin'] = pd.cut(df.mc_energy, bins)

    x = df.mc_energy.values
    y = distance

    bin_means, bin_edges, binnumber = binned_statistic(x, y, statistic='median', bins=bins)

    b_68, bin_edges, binnumber = binned_statistic(x, y, statistic=lambda y: np.percentile(y, 68), bins=bins)
    b_32, bin_edges, binnumber = binned_statistic(x, y, statistic=lambda y: np.percentile(y, 32), bins=bins)

    bin_centers = np.sqrt(bin_edges[1:] * bin_edges[:-1])

    plt.scatter(x, y, s=0.2, color='#909090')
    plt.hlines(bin_means, bin_edges[:-1], bin_edges[1:], lw=1, colors=color)
    c = bin_centers
    c[0] = bin_edges.min()
    c[-1] = bin_edges.max()
    plt.fill_between(c, b_32, b_68, facecolor=color, alpha=0.4, interpolate=False)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('Distance to True Position / degree')
    plt.xlabel('Energy / TeV')
    plt.ylim([0.001, 100.8])
    plt.xlim([bin_edges.min(), bin_edges.max()])
    if title:
        plt.title(title)
    plt.tight_layout()
    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == "__main__":
    main()
