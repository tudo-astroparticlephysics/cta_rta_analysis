import click
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from sklearn.metrics import roc_auc_score
from tqdm import tqdm
import fact.io

columns = ['array_event_id', 'gamma_prediction', 'telescope_type_name']


def add_rectangles(ax, offset=0.1):

    kwargs = {
        'linewidth': 1,
        'edgecolor': 'white',
        'facecolor': 'white',
        'alpha': 0.6,
    }

    # left
    rect = patches.Rectangle((0 - offset, 0), 0 + offset, 1 + offset, **kwargs)
    ax.add_patch(rect)

    # right
    rect = patches.Rectangle((1, 0 - offset), 0 + offset, 1 + offset, **kwargs)
    ax.add_patch(rect)

    # top
    rect = patches.Rectangle((0, 1), 1 + offset, 0 + offset, **kwargs)
    ax.add_patch(rect)

    # bottom
    rect = patches.Rectangle((0 - offset, 0), 1 + offset, 0 - offset, **kwargs)
    ax.add_patch(rect)



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
    '-o', '--output', type=click.Path(
        exists=False,
        dir_okay=False,
    ))
def main(predicted_gammas, predicted_protons, output):
    telecope_events = fact.io.read_data(predicted_gammas, key='telescope_events', columns=columns).dropna()
    array_events = fact.io.read_data(predicted_gammas, key='array_events', columns=['array_event_id', 'mc_energy'])
    gammas = pd.merge(telecope_events, array_events, on='array_event_id')

    telecope_events = fact.io.read_data(predicted_protons, key='telescope_events', columns=columns).dropna()
    array_events = fact.io.read_data(predicted_protons, key='array_events', columns=['array_event_id', 'mc_energy'])
    protons = pd.merge(telecope_events, array_events, on='array_event_id')


    bin_edges = np.logspace(-2, 2, 50)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_widths = np.diff(bin_edges)

    gammas['energy_bin'] = pd.cut(gammas.mc_energy, bin_edges)
    protons['energy_bin'] = pd.cut(protons.mc_energy, bin_edges)

    for tel_type in ['SST', 'MST', 'LST']:
        aucs = []
        for b in tqdm(gammas.energy_bin.cat.categories):

            tel_gammas = gammas[(gammas.energy_bin == b) & (gammas.telescope_type_name == tel_type)]
            tel_protons = protons[(protons.energy_bin == b) & (protons.telescope_type_name == tel_type)]

            if len(tel_gammas) < 30 or len(tel_protons) < 30:
                aucs.append(0)
            else:
                mean_prediction_gammas = tel_gammas.groupby('array_event_id')['gamma_prediction'].mean()
                gamma_labels = np.ones_like(mean_prediction_gammas)

                mean_prediction_protons = tel_protons.groupby('array_event_id')['gamma_prediction'].mean()
                proton_labels = np.zeros_like(mean_prediction_protons)

                y_score = np.hstack([mean_prediction_gammas, mean_prediction_protons])
                y_true = np.hstack([gamma_labels, proton_labels])

                aucs.append(roc_auc_score(y_true, y_score))


        plt.errorbar(
            bin_centers,
            aucs,
            xerr=bin_widths / 2.0,
            linestyle='',
            label=tel_type
        )

    plt.ylim([0, 1])
    plt.xscale('log')
    plt.legend()
    # add_rectangles(plt.gca())

    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == '__main__':
    main()