from dask_drmaa import DRMAACluster
import click
from dask.distributed import Client
from process_simtel import process_file, verify_file
import fact.io


@click.command()
@click.argument(
    'input_files', type=click.Path(
        exists=True,
        dir_okay=False,
    ), nargs=-1)
@click.argument(
    'output_file', type=click.Path(
        dir_okay=False,
    ))
@click.option('-n', '--n_workers', default=5)
def main(input_files, output_file, n_workers):
    cluster = DRMAACluster()
    client = Client(cluster)
    cluster.start_workers(n_workers)

    f = client.map(process_file, input_files)
    frames = client.gather(f)
    runs = [r[0] for r in frames]
    array_events = [r[1] for r in frames]
    telescope_events = [r[2] for r in frames]

    fact.io.write_data(runs, output_file, key='runs', mode='a')
    fact.io.write_data(array_events, output_file, key='array_events', mode='a')
    fact.io.write_data(telescope_events, output_file, key='telescope_events', mode='a')

    verify_file(output_file)


if __name__ == '__main__':
    main()
