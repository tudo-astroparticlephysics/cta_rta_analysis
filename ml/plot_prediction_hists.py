import click
import numpy as np
import matplotlib.pyplot as plt
import fact.io
from cycler import cycler


@click.command()
@click.argument(
    'predicted_gammas', type=click.Path(
        exists=True,
        dir_okay=False,
    ))
@click.argument(
    'predicted_protons', type=click.Path(
        exists=True,
        dir_okay=False,
    ))
@click.option(
    '-o', '--output_file', type=click.Path(
        exists=False,
        dir_okay=False,
    ))
def main(predicted_gammas, predicted_protons, output_file):
    bins = np.linspace(0, 1, 100)

    gammas = fact.io.read_data(predicted_gammas, key='telescope_events')

    fig, ax = plt.subplots(1)
    for name, group in gammas.groupby('telescope_type_name'):
        ax.hist(group.gamma_prediction.values, bins=bins, label=f'gamma prediction {name}', histtype='step', linewidth=2)


    color_cycle = cycler(color=['gray', 'darkgray', 'black'])
    ax.set_prop_cycle(color_cycle)
    protons = fact.io.read_data(predicted_protons, key='telescope_events')
    for name, group in protons.groupby('telescope_type_name'):
        ax.hist(group.gamma_prediction.values, bins=bins, label=f'proton prediction {name}', histtype='step')


    plt.legend(loc='upper left')

    if output_file:
        plt.savefig(output_file)
    else:
        plt.show()


if __name__ == '__main__':
    main()
