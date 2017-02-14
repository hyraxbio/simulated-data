from seq2simulate import run_make_proviral_data

import click
import os
import appdirs
import uuid

def get_working_dir():
    """
    Returns a default working directory.
    """

    return os.path.join(
        appdirs.user_data_dir('tests', 'simulated'), 
        str(uuid.uuid4())
    )

@click.command()

@click.option(
    '--sequences', 
    help='Path to a FASTA file of sequence(s).', 
    required=True,
    type=click.Path(exists=True, readable=True, resolve_path=True)
)

@click.option(
    '--out', 
    help='Path to the desired output folder (will be created).',
    type=click.Path(resolve_path=True)
)

@click.option(
    '--platform',
    help='Sequencing platform to simulate.',
    required=True,
    type=click.Choice(['illumina', 'ion', 'pacbio-clr', 'pacbio-ccs', 'roche'])
)

@click.option(
    '--paired-end',
    help='Paired-end sequencing.  Only valid with illumina data.',
    is_flag=True
)

@click.option(
    '--working-dir',
    help='Returns a specific working directory location. ' \
         'Defaults to a randomly generated uuid in an appropriate ' \
         'data folder.', 
    type=click.Path(resolve_path=True)
)

@click.option(
    '--proviral-fraction', 
    help='Fraction of final population to be hypermutated.', 
    default=0.1,
)

def run(
    sequences, out, platform, working_dir, paired_end, proviral_fraction,
):
    """
    """
    if paired_end and platform != "illumina":
        raise ValueError(
            'Only illumina paired-end reads can be simulated.'
        )

    if working_dir is None:
        working_dir = get_working_dir()

    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    if out is None:
        out = working_dir
    else:
        if not os.path.isdir(out):
            os.makedirs(out)

    if not 0 <= proviral_fraction <= 1:
        raise ValueError('proviral_fractions must be in range 0-1')

    run_make_proviral_data.run_proviral(sequences, working_dir, out, platform, paired_end, proviral_fraction)

if __name__ == '__main__':
    run()