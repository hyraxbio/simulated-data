import random
import sys
from Bio import SeqIO

import sequence



def parse_sequences(raw_susceptible_file, raw_resistant_file):
    """
    Parses a raw sequence file into a Sequence object.

    Args:
        raw_sequence_file: the filename where sequences are stored.
        raw_drm_file: a file containing a comma-separated list of DRM 
        positions.

    Returns:
        A list of Sequence objects.
    """

    print 'Parsing sequences.'

    sequences = []
    for susceptible, resistant in \
        zip(SeqIO.parse(raw_susceptible_file, "fasta"), \
            SeqIO.parse(raw_resistant_file, "fasta")):
        
        sequences.append(sequence.Sequence(
                            susceptible, 
                            resistant
                        ))
        print '.',
        sys.stdout.flush()

    print
    return sequences

            

