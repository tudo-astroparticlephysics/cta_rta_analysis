import click
import matplotlib.pyplot as plt
import fact.io
import numpy as np


@click.command()
@click.argument('diffuse', type=click.Path(exists=True))
@click.argument('point', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
@click.option('-t', '--title', default=None)
def main(diffuse, point, output, title):
    fig, [left, right] = plt.subplots(2, 2, figsize=(7, 10), constrained_layout=True)
    for l, file_path, axs in zip(['diffuse gammas', 'point-like gammas'], [diffuse, point], [left, right]):
        df = fact.io.read_data(file_path, key='array_events').dropna()
        x = df.core_x_prediction
        y = df.core_y_prediction

        mc_x = df.mc_core_x
        mc_y = df.mc_core_y

        ax1 = axs[0]
        bins = np.linspace(-700, 700, 40)
        ax1.hist2d(mc_x, mc_y, bins=bins, cmap='gray')
        ax1.set_xlabel('True X')
        ax1.set_ylabel('True Y')
        ax1.set_aspect('equal')

        ax2 = axs[1]
        ax2.hist2d(x, y, bins=bins, cmap='gray')
        ax2.set_aspect('equal')
        ax2.set_xlabel('Reconstructed X')
        ax2.set_ylabel('Reconstructed Y')

    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == "__main__":
    main()
