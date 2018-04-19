import click
import matplotlib.pyplot as plt
import numpy as np
from sensitivity import read_sensitivity_fits



@click.command()
@click.argument('input_fits', type=click.Path(exists=True), nargs=-1)
@click.option('-l', '--label', multiple=True)
@click.option('-c', '--color', multiple=True)
@click.option('-o', '--output', type=click.Path(exists=False))
def main(
    input_fits, label, color,
    output,
):
    '''
    Plots cuts used in a sensitivity curve.
    '''


    if label and len(input_fits) != len(label):
        print('Must pass as many labels as gamma files as proton files')

    if color and len(label) != len(color):
        print('Must pass as many colors as labels')


    if not label:
        label = [None] * len(input_fits)

    if not color:
        color = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 6))
    for fits_file, l, c in zip(input_fits, label, color):
        table, bin_edges, bin_center, bin_widths = read_sensitivity_fits(fits_file)

        x, y = bin_center.value, np.sqrt(table['radius'].value)
        xerr = bin_widths.value * 0.5

        y = np.ma.masked_where(y <= 0.05, y)

        ax1.errorbar(
            x, y,
            xerr=xerr,
            marker='.',
            linestyle='--',
            capsize=0,
            ecolor='gray',
            color=c,
            ms=0,
        )

        y = table['threshold']
        y = np.ma.masked_where(y <= 0, y)
        ax2.errorbar(
            x, y,
            xerr=xerr,
            marker='.',
            linestyle='--',
            ecolor='gray',
            capsize=0,
            label=l,
            ms=0,
            color=c,
        )

    ax1.set_ylabel('on region radius / degree')
    ax1.set_xlabel('Energy / TeV')
    ax1.set_xscale('log')
    ax1.set_ylim([0, 0.5])

    ax2.set_ylabel('prediction threshold')
    ax2.set_xlabel('Energy / TeV')
    ax2.set_xscale('log')
    ax2.set_ylim([0, 1])


    if label[0]:
        plt.legend()

    if output:
        plt.savefig(output)
    else:
        plt.show()



if __name__ == '__main__':
    main()
