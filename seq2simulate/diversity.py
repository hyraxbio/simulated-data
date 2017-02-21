import os
import random
import sys
import shutil
import time

import Bio

import evolveagene
import sierra_ws as sierra
from hypermutation import hypermutate
import hiv_drms

output_filename = 'evolveagene_checked.fasta'

hiv_pol_dnds = 0.23
hiv_pol_substitution_rate = 0.00072
hiv_pol_lambda = 0.00001 #probability of indel at locus (~codon) in branch
hiv_pol_ti_td = 0.1 #ratio of insertions to deletions
min_taxa_to_keep = 4
max_tries = 3

# proviral hypermutations per hundred bp
hypermutation_rate = 3

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


def simulate(sequence, working_dir, num_taxa=10, hypermutate_seqs=False, include_deletions=False):
    """
    Produce a simulated set of sequences that contain no added or removed DRMs
    with respect to the original sequence. Also, optionally hypermutate evolved
    sequences.

    Args:
        sequence: The sequence to simulate.
        working_dir: Where the simulated sequence files go.
        num_taxa: The number of taxa to evolve.
        hypermutate_seqs: Perform hypermutation on the sequence.

    Returns:
        A tuple of filenames of the susceptible and resistant files.
    """
    

    sequences = {
        ('susceptible', sequence.susceptible),
        ('resistant', sequence.resistant)
    }

    filenames = {}

    for name, seq in sequences:
        allowed_sequences = _simulate_evolution(name, seq, sequence, working_dir, num_taxa)
        if hypermutate_seqs:
            allowed_sequences = _simulate_hypermutation(allowed_sequences)
        olds = [len(i.seq) for i in allowed_sequences]
        if include_deletions:
            allowed_sequences = _simulate_deletions(allowed_sequences, freq=0.4)
        news = [len(i.seq) for i in allowed_sequences]
        diffs = [i!=j for i,j in zip(olds, news)]
        print('freq:', float(sum(diffs))/len(diffs))

        full_filename = os.path.join(
            working_dir, 
            seq.id + "_" + name + "_" + output_filename
        )

        Bio.SeqIO.write(allowed_sequences, full_filename, 'fasta')
        filenames[name] = full_filename

    return filenames

def _simulate_evolution(name, seq, sequence, working_dir, num_taxa):

    try_idx = 0
    while True:
        branch_length = (
            hiv_pol_substitution_rate \
                * sequence.years_infected \
                * random.uniform(0.6, 1.0)
        )

        evolve_file = evolveagene.run( seq, num_taxa, branch_length,
            hiv_pol_dnds, hiv_pol_lambda, hiv_pol_ti_td,
            working_dir)

        evolved_sequences = [s for s in Bio.SeqIO.parse(evolve_file, 'fasta')]

        drms = []
        if name == 'resistant':
            drms = sequence.drms

        print "Checking taxa for missing or introduced DRMs."
        print "Years infected:", sequence.years_infected

        try:
            allowed_sequences = [
                s for s in evolved_sequences \
                    if drms_unchanged(seq.id, drms, 
                        sierra.get_drms(s, known_drms=hiv_drms))
            ]
        except:
            allowed_sequences = []

        if len(allowed_sequences) >= min_taxa_to_keep:

            print "Kept a total of", len(allowed_sequences), "evolved " \
                "sequences."
            
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

    return allowed_sequences


def _simulate_hypermutation(sequences):

    print('\n-------------------------------------------------------')
    print('Hypermutating evolved sequences (rate = {} per 100 bp).'.format(hypermutation_rate))
    print('-------------------------------------------------------\n')
    hypermutation_rates = [int(hypermutation_rate * len(str(s.seq)) // 100) for s in sequences]
    hyper_evolved_sequences = hypermutate.mutate_sequences([str(s.seq) for s in sequences], hypermutation_rates)
    for i, hseq in enumerate(hyper_evolved_sequences):
        sequences[i].seq = Bio.Seq.Seq(hseq, alphabet=Bio.Alphabet.SingleLetterAlphabet())
    return sequences

def _simulate_deletions(sequences, freq=0.4, strip_deletions=True):
    """
    Args:
        sequences: list of DNA strings
        freq: probability of deletion

    Returns:
        list of sequences

    Bruner (2016) proviral defects:
        2%  - intact
        5%  - packaging signal and major splice donor site deletion (15-97 bp)
        20% - deletion within 5' half of genome (15-6000 bp) 
        35% - deletion within 3' half of genome (15-6000 bp)
        20% - very large internal deletion (6000-8000 bp)
        7%  - hypermutated
        8%  - hypermutated and deleted
        3%  - insertions and inversions

    Assuming that large deletions and packaging signal/major splice donor site
    deletions don't apply to individual gene sequences:

        10%  - intact
        35%  - hypermutated
        40%  - hypermutated and deleted
        15%  - insertions and inversions

    Refs:        
        Bruner, K. M., Murray, A. J., Pollack, R. A., Soliman, M. G., Laskey, S.
        B., Capoferri, A. A., … Siliciano, R. F. (2016). Defective proviruses rapidly
        accumulate during acute HIV-1 infection. Nature Medicine, 22(August),
        doi:10.1038/nm.4156. https://doi.org/10.1038/nm.4156

    """
    print('\n-------------------------------------------------------')
    print('Making deletions in evolved sequences (probability = {}).'.format(freq))
    print('-------------------------------------------------------\n')

    DEL_SIGN = '-'
    if strip_deletions:
        DEL_SIGN = ''

    string_seqs = [str(s.seq) for s in sequences]
    for i, seq in enumerate(string_seqs):
        if random.uniform(0, 1) <= freq:
            del_ranges = sorted([random.randrange(0, len(seq)), random.randrange(0, len(seq))])
            string_seqs[i] = seq[0:del_ranges[0]] + DEL_SIGN*(del_ranges[1]-del_ranges[0]) + seq[del_ranges[0]:del_ranges[1]]
    for i, sseq in enumerate(sequences):
        sequences[i].seq = Bio.Seq.Seq(string_seqs[i], alphabet=Bio.Alphabet.SingleLetterAlphabet())
    return sequences
