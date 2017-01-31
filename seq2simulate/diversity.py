import os
import random
import sys
import shutil
import time

import Bio

import evolveagene
import sierra_wrapper as sierra

output_filename = 'evolveagene_checked.fasta'

hiv_pol_dnds = 0.23
hiv_pol_substitution_rate = 0.00072
tree_type = evolveagene.RANDOM
selection_type = evolveagene.CONSTANT
num_taxa = 10
min_taxa_to_keep = 4
max_tries = 3


def drms_unchanged(id, drms1, drms2):
    """
    Compare two DRM lists for equality.

    Args:
        drms1, drms2: The two lists to compare.

    Returns:
        True if equal, else False.
    """

    if sorted(drms1) != sorted(drms2):
        print drms1, "!=", drms2
    return sorted(drms1) == sorted(drms2)


def simulate(sequence, working_dir):
    """
    Produce a simulated set of sequences that contain no added or removed
    DRMs with respect to the original sequence.

    Args:
        sequence: The sequence to simulate.
        working_dir: Where the simulated sequence files go.

    Returns:
        A tuple of filenames of the susceptible and resistant files.
    """

    sequences = {
        ('susceptible', sequence.susceptible),
        ('resistant', sequence.resistant)
    }

    filenames = {}

    for name, seq in sequences:
        out_file = os.path.join(working_dir, output_filename)

        try_idx = 0
        while True:
            branch_length = (
                hiv_pol_substitution_rate \
                    * sequence.years_infected \
                    * random.uniform(0.6, 1.0)
            )

            evolve_file = evolveagene.run(
                seq, tree_type, selection_type, num_taxa, 
                branch_length, hiv_pol_dnds, working_dir
            )

            ## evolveagene OOMs aggressively, which is a pain
            #tries = 10
            #while not os.path.isfile(evolve_file) and tries > 0:
            #    # doze a bit so the threads get out of sync
            #    time.sleep(random.randint(1, 5))
            #    evolve_file = evolveagene.run(
            #        seq, tree_type, selection_type, num_taxa, 
            #        branch_length, hiv_pol_dnds, working_dir
            #    )
            #    tries -= 1

            drms = []
            if name == 'resistant':
                drms = sequence.drms

            print "Checking taxa for missing or introduced DRMs."
            print "Years infected:", sequence.years_infected

            try:
                allowed_sequences = [
                    s for s in Bio.SeqIO.parse(evolve_file, 'fasta') \
                        if drms_unchanged(seq.id, drms, 
                            sierra.get_drms(s))
                ]
            except:
                allowed_sequences = []

            if len(allowed_sequences) >= min_taxa_to_keep:

                print "Kept a total of", len(allowed_sequences), "evolved " \
                    "sequences."
                
                full_filename = os.path.join(
                    working_dir, 
                    seq.id + "_" + name + "_" + output_filename
                )

                Bio.SeqIO.write(allowed_sequences, full_filename, 'fasta')
                filenames[name] = full_filename

                os.unlink(evolve_file)
                break
            elif try_idx == max_tries - 1:
                sequence.years_infected -= 1
                sequence.years_infected = max(1, sequence.years_infected)
                try_idx = 0
                print "Lowering year count, evolving again."
            else:
                try_idx += 1
                print "Too many sequences have DRMs, evolving again."

    return filenames
