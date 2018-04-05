import click
import fact.io
from colorama import Fore


@click.command()
@click.argument(
    'input_file_path', type=click.Path(
        exists=True,
        dir_okay=False,
    ))
def main(input_file_path):
    runs = fact.io.read_data(input_file_path, key='runs')
    runs.set_index('run_id', drop=True, verify_integrity=True, inplace=True)

    telescope_events = fact.io.read_data(input_file_path, key='telescope_events')
    telescope_events.set_index('telescope_event_id', drop=True, verify_integrity=True, inplace=True)

    array_events = fact.io.read_data(input_file_path, key='array_events')
    array_events.set_index('array_event_id', drop=True, verify_integrity=True, inplace=True)

    assert len(array_events) == len(telescope_events.array_event_id.unique())

    assert len(runs) == len(telescope_events.run_id.unique())
    assert len(runs) == len(array_events.run_id.unique())
    print(Fore.GREEN + 'File looks OK.')


if __name__ == '__main__':
    main()
