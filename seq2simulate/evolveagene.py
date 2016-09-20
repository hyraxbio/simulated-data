import os
import sys
import uuid
from Bio import SeqIO
from subprocess import call

fasta_output_filename = 'evolveagene_raw'
evolveagene_extension = '_Unaligned.FASTA'

CONSTANT, PURIFYING, POSITIVE, BOTH = range(4)
SYMMETRIC, RANDOM = range(2)

tree_types = {
    SYMMETRIC: 'sym',
    RANDOM: 'ran'
}

selection_types = {
    CONSTANT: 'con',
    PURIFYING: 'pur',
    POSITIVE: 'pos',
    BOTH: 'both'
}

def run(
    sequence, tree_type, selection_type, 
    num_taxa, branch_length, dnds, working_dir):
    """
    Run the EvolveAGene4 sequence simulator.

    Args:
        sequence: A BioPython sequence to simulate.
        tree_type: Either SYMMETRIC or RANDOM
        selection_type: CONSTANT, PURIFYING, POSITIVE or BOTH.
        num_taxa: The number of sequences to simulate.
        branch_length: The average substitution rate of the branches.
        dnds: The chance that a nonsynonymous mutation will be accepted.
        working_dir: The director to which to save the generated file.

    Returns:
        True on completion.
    """

    print 'Evolving sequence', sequence.id, '...',
    sys.stdout.flush()

    sequence_file = os.path.join(working_dir, 
        fasta_output_filename + str(uuid.uuid4()))

    # EvolveAGene doesn't like fasta headers, so we write plain text
    # sequence rather than using SeqIO
    with open(sequence_file, 'w') as out:
        out.write(str(sequence.seq))

    arg_list = ['EvolveAGene4',
        '-f', '\"%s\"' % sequence_file, 
        '-t', tree_types[tree_type],
        '-n', str(num_taxa),
        '-a', str(dnds),
        '-b', str(branch_length),
        '-ss', selection_types[selection_type],
        ]

    call(" ".join(arg_list), shell=True)

    print 'evolved.'


    return sequence_file + evolveagene_extension
