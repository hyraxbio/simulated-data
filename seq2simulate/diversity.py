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
# hypermutation_rate = 3

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


def simulate(sequence, working_dir, num_taxa=10, 
             hypermutate_seqs=False, 
             include_deletions=False,
             include_insertions=False,
             include_frameshifts=False,
             include_stop_codons=False,
             include_inversions=False,
            ):
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
        allowed_sequences_strings, _ = _convert_seqs_to_strs(allowed_sequences)
        mutation_types  = [hypermutate_seqs, 
                           include_deletions, 
                           include_insertions, 
                           include_frameshifts, 
                           include_stop_codons, 
                           include_inversions]
        mutation_funcs = [_simulate_hypermutation, 
                          _simulate_deletions, 
                          _simulate_insertions, 
                          _simulate_frameshifts, 
                          _simulate_stop_codons, 
                          _simulate_inversions]
        for mutation_type, mutation_func in zip(mutation_types, mutation_funcs):
            if mutation_type:
                allowed_sequences_strings, seq_diffs = mutation_func(allowed_sequences_strings)
        _update_seq_reads_from_strs(allowed_sequences, allowed_sequences_strings)

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


def _simulate_hypermutation(sequences, hypermutation_rate=3):

    print('\n-------------------------------------------------------')
    print('Hypermutating evolved sequences (rate = {} per 100 bp).'.format(hypermutation_rate))
    print('-------------------------------------------------------\n')
    hypermutation_rates = [int(hypermutation_rate * len(s) // 100) for s in sequences]
    hyper_evolved_sequences = hypermutate.mutate_sequences(sequences, hypermutation_rates)
    seq_diffs = [get_diffs1(sequences[i], hyper_evolved_sequences[i]) for i in range(len(sequences))]

    return sequences, seq_diffs

def _simulate_deletions(sequences, freq=0.4, strip_deletions=True, max_length=120, min_length=15, no_frameshifts=True):
    """
    Args:
        sequences: list of DNA strings
        freq: probability of deletion
        strip_deletions: do not include '-' characters for deleted residues
        max_length: maximum length of deletion
        min_length: minimum length of deletion
        no_frameshifts: do not introduce frameshifts

    Returns:
        list of sequences

    """
    print('\n-------------------------------------------------------')
    print('Making deletions in evolved sequences (probability = {}).'.format(freq))
    print('-------------------------------------------------------\n')

    DEL_SIGN = '-'
    if strip_deletions:
        DEL_SIGN = ''

    seq_diffs = []
    for i, seq in enumerate(sequences):
        if min_length > len(seq):
            min_length = len(seq)//4
        if max_length > len(seq):
            max_length = len(seq)//2
        if random.uniform(0, 1) <= freq:
            if no_frameshifts: 
                del_start = random.randrange(0, (len(seq)-min_length)//3)*3
                del_length = random.randint(0, min(max_length, len(seq)-del_start)//3)*3
            else:
                del_start = random.randrange(0, len(seq)-min_length)
                del_length = random.randint(0, min(max_length, len(seq)-del_start))
            sequences[i] = seq[0:del_start] + DEL_SIGN*del_length + seq[del_start+del_length:]
            seq_diffs.append([del_start, del_start+del_length])
        else:
            seq_diffs.append([])
    return sequences, seq_diffs

def _simulate_insertions(sequences, freq=0.2, max_length=100, min_length=15, no_frameshifts=True):
    """
    Args:
        sequences: list of DNA strings
        freq: probability of insertion
        max_length: of random insertion 
        min_length: of random insertion
        no_frameshifts: do not introduce frameshifts

    Returns:
        list of sequences


    """
    print('\n---------------------------------------------------------')
    print('Making insertions in evolved sequences (probability = {}).'.format(freq))
    print('---------------------------------------------------------\n')

    seq_diffs = []
    for i, seq in enumerate(sequences):
        if min_length > len(seq):
            min_length = len(seq)//4
        if max_length > len(seq):
            max_length = len(seq)//2
        if random.uniform(0, 1) <= freq:
            if no_frameshifts:
                ins_start = random.randrange(0, (len(seq)-min_length)//3)*3
                ins_length = random.randint(0, min(max_length, len(seq)-ins_start)//3)*3
            else:
                ins_start = random.randrange(0, len(seq)-min_length)
                ins_length = random.randint(0, min(max_length, len(seq)-ins_start))
            sequences[i] = seq[0:ins_start] + ''.join([random.choice('ATGC') for insertion in range(ins_length)]) + seq[ins_start:]
            seq_diffs.append([ins_start, ins_start+ins_length])
        else:
            seq_diffs.append([])
    return sequences, seq_diffs

def _simulate_frameshifts(sequences, freq=0.3, strip_deletions=True):
    """
    Args:
        sequences: list of DNA strings
        freq: probability of deletion

    Returns:
        list of sequences

    """
    print('\n-------------------------------------------------------')
    print('Making frameshifts in evolved sequences (probability = {}).'.format(freq))
    print('-------------------------------------------------------\n')

    DEL_SIGN = '-'
    if strip_deletions:
        DEL_SIGN = ''

    seq_diffs = []
    for i, seq in enumerate(sequences):
        if random.uniform(0, 1) <= freq:
            del_start = random.randrange(0, len(seq))
            sequences[i] = seq[0:del_start] + DEL_SIGN + seq[del_start+1:]
            seq_diffs.append([del_start])
        else:
            seq_diffs.append([])
    return sequences, seq_diffs

def _simulate_stop_codons(sequences, freq=0.5):
    """
    Args:
        sequences: list of DNA strings
        freq: probability of mutation

    Returns:
        list of sequences

    """
    print('\n-------------------------------------------------------')
    print('Making stop codons in evolved sequences (probability = {}).'.format(freq))
    print('-------------------------------------------------------\n')

    stop_codons = ['TAG', 'TAA', 'TGA']
    seq_diffs = []
    for i, seq in enumerate(sequences):
        if random.uniform(0, 1) <= freq:
            codon_inds = [k*3 for k in range(len(seq)//3)]
            putative_codons = []
            for ind in codon_inds:
                codon = seq[ind:ind+3]
                stop_codon, score = _closest_match(codon, stop_codons)
                if score == 2:
                    putative_codons.append([ind, stop_codon]) 
            if len(putative_codons) > 0:
                codon_replacement = random.choice(putative_codons)
                sequences[i] = seq[0:codon_replacement[0]] + codon_replacement[1] + seq[codon_replacement[0]+3:]
                seq_diffs.append([codon_replacement[0]])
        else:
            seq_diffs.append([])
    return sequences, seq_diffs

def _closest_match(s1, s2):
    scores = [_sim_score(s1, i) for i in s2]
    return s2[scores.index(max(scores))], max(scores)

def _sim_score(s1, s2):
    return sum([i==j for i,j in zip(s1, s2)])

def _simulate_inversions(sequences, freq=0.3, max_length=100, min_length=5):
    """
    Args:
        sequences: list of DNA strings
        freq: probability of insertion
        max_length: of random insertion 
        min_length: of random insertion

    Returns:
        list of sequences


    """
    print('\n---------------------------------------------------------')
    print('Making inversions in evolved sequences (probability = {}).'.format(freq))
    print('---------------------------------------------------------\n')

    seq_diffs = []
    for i, seq in enumerate(sequences):
        if random.uniform(0, 1) <= freq:
            ins_start = random.randrange(0, len(seq)-min_length)
            ins_length = random.randint(0, min(max_length, len(seq)-ins_start))
            sequences[i] = seq[0:ins_start] + seq[ins_start:ins_start + ins_length][::-1] + seq[ins_start + ins_length:]
            seq_diffs.append([ins_start, ins_start + ins_length])
        else:
            seq_diffs.append([])
    return sequences, seq_diffs

def _convert_seqs_to_strs(sequences):
    return [str(s.seq) for s in sequences], [s.id for s in sequences]

def _update_seq_reads_from_strs(sequences, string_seqs):
    for i, sseq in enumerate(sequences):
        sequences[i].seq = Bio.Seq.Seq(string_seqs[i], alphabet=Bio.Alphabet.SingleLetterAlphabet())

def get_diffs1(seq0, seq1):
    """
    Simply return a list of 1-based different indices between two strings.
    """
    assert len(seq0) == len(seq1)
    return [i for i in range(len(seq0)) if seq0[i] != seq1[i]]
