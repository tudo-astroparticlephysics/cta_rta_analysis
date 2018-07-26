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


    df = fact.io.read_data(diffuse, key='array_events').dropna()
    x = df.core_x_prediction
    y = df.core_y_prediction

    mc_x = df.mc_core_x
    mc_y = df.mc_core_y

    d = (mc_x - x)**2 + (mc_y - y)**2
    bins = np.linspace(0, 100000, 40)
    plt.hist(d, bins=bins, histtype='step', density=True, lw=3)

    df = fact.io.read_data(point, key='array_events').dropna()
    x = df.core_x_prediction
    y = df.core_y_prediction

    mc_x = df.mc_core_x
    mc_y = df.mc_core_y

    d = (mc_x - x)**2 + (mc_y - y)**2
    plt.hist(d, bins=bins, histtype='step', density=True, lw=3)



    plt.xlabel('distance to true impact squared')



    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == "__main__":
    main()
