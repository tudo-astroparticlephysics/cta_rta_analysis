import click
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import cartesian_to_spherical, Angle, SkyCoord, EarthLocation
from dateutil import parser
import seaborn as sns
from scipy.stats import binned_statistic
from matplotlib.colors import LogNorm


@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_file', type=click.Path(exists=False))
def main(input_file, output_file):
    df = pd.read_csv(input_file).dropna()

    x = df['stereo:estimated_direction:x']
    y = df['stereo:estimated_direction:y']
    z = df['stereo:estimated_direction:z']

    r, lat, lon = cartesian_to_spherical(x.values * u.m, y.values * u.m, z.values * u.m)

    alt = Angle(90 * u.deg - lat)
    mc_alt = Angle(df['mc:alt'].values * u.rad)

    az = Angle(lon).wrap_at(180 * u.deg)
    mc_az = Angle(df['mc:az'].values * u.rad).wrap_at(180 * u.deg)

    paranal = EarthLocation.of_site('paranal')
    dt = parser.parse('1987-09-20 22:15')

    c = SkyCoord(
        alt=alt,
        az=az,
        obstime=dt,
        frame='altaz',
        location=paranal,
    )

    c_mc = SkyCoord(
        alt=mc_alt,
        az=mc_az,
        obstime=dt,
        frame='altaz',
        location=paranal,
    )

    df['spherical_distance'] = c.separation(c_mc)

    df['mc:energy'] = df['mc:energy'].apply(lambda x: np.log10(x * 1000))
    #

    df = df[df['stereo:estimated_direction:z'] > 0]

    bins = np.linspace(df['mc:energy'].min(), df['mc:energy'].max(), 10)
    df['energy_bin'] = pd.cut(df['mc:energy'], bins)
    #
    x = df['mc:energy'].values
    #
    y = df['spherical_distance'].values
    #

    bin_means, bin_edges, binnumber = binned_statistic(x, y, statistic='median', bins=25)
    b_68, bin_edges, binnumber = binned_statistic(x, y, statistic=lambda y: np.percentile(y, 68), bins=25)
    b_32, bin_edges, binnumber = binned_statistic(x, y, statistic=lambda y: np.percentile(y, 32), bins=25)

    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = (bin_edges[1:] + bin_edges[0:-1])/2

    # plt.errorbar(bin_centers, bin_means, yerr=[b_95, b_05], xerr=bin_width/2, elinewidth=2, lw=0, ecolor='#348ABD')
    plt.scatter(x, y, s=0.2, color='#909090')
    plt.hlines(bin_means, bin_edges[:-1], bin_edges[1:], lw=1, colors='#0a4f88')
    c = bin_centers
    c[0] = bin_edges.min()
    c[-1] = bin_edges.max()
    plt.fill_between(c, b_32, b_68, facecolor='#0a4f88', alpha=0.4,  interpolate=True)
    # plt.hist2d(df['mc:energy'], df['spherical_distance'], bins=100, norm=LogNorm())
    # whis=[15.7, 15.7+68.2]
    # sns.boxplot(x='energy_bin', y='spherical_distance', data=df, fliersize=0.5, linewidth=1)
    plt.yscale('log')
    plt.ylabel('Distance to True Position / degree')
    plt.xlabel('log(Energy / TeV)')
    plt.ylim([0, 0.8])
    plt.xlim([bin_edges.min(), bin_edges.max()])

    plt.tight_layout()
    plt.savefig(output_file)


if __name__ == "__main__":
    main()
