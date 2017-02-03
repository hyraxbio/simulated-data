import os
import sys
import uuid
from Bio import SeqIO
from subprocess import call

from revolver.evolve import evolve

fasta_output_filename = 'evolveagene_raw'

def run(
    sequence, num_taxa, branch_length, dnds, lmbda, ti_td, working_dir):
    """
    Run the EvolveAGene4 sequence simulator.

    Args:
        sequence: A BioPython sequence to simulate.
        num_taxa: The number of sequences to simulate.
        branch_length: The average substitution rate of the branches (by
        default this is number of nucleotide substitutions per codon)
        dnds: Ratio of non-synonymous to synonymous mutations
        lmbda: probability of indel at locus/codon
        ti_td: ratio of insertions to deletions
        working_dir: The director to which to save the generated file.

    Returns:
        True on completion.
    """

    print 'Evolving sequence', sequence.id, '...',
    sys.stdout.flush()

    new_sequences = evolve(str(sequence.seq), 
        taxa=num_taxa, 
        t=branch_length,
        omega=dnds, 
        lmbda=lmbda,
        ti_td=ti_td,
        codon_freq='F1x4',
        log=False,
        scale_q=True,
        strip_deletions=False,
        verbose=True)  


    sequence_file = os.path.join(working_dir, 
        fasta_output_filename + str(uuid.uuid4()))

    print 'evolved.'
    with open(sequence_file, 'w') as out:
        for i,j in enumerate(new_sequences):
            out.write('>{}\n'.format(i))
            out.write(j+'\n')
    
    return sequence_file
