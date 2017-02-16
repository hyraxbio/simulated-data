from Bio import SeqIO, Seq, Alphabet
import os
import random
import json

from hypermutation import hypermutate
import platform as plat
import art
from custom_generators import parse_fastq

# proviral hypermutations per hundred bp
hypermutation_rate = 3

def get_diffs(seq0, seq1):
    """
    Simply return a list of different indices between two strings.
    """
    assert len(seq0) == len(seq1)
    return [i for i in range(len(seq0)) if seq0[i] != seq1[i]]

def hypermutate_sequences(sequences, working_dir):
    """
    Hypermutates a set of sequence strings.

    Args:
        sequences: list of DNA strings
        working_dir: temporary directory for storing output files 

    Returns:
        FASTA filename of hypermutated sequences
    """
    print('\n-------------------------------------------------------')
    print('Hypermutating evolved sequences (rate = {} per 100 bp).'.format(hypermutation_rate))
    print('-------------------------------------------------------\n')
    string_seqs = [str(s.seq) for s in sequences]
    hypermutation_rates = [int(hypermutation_rate * len(str(s.seq)) // 100) for s in sequences]
    hyper_evolved_sequences = hypermutate.mutate_sequences(string_seqs, hypermutation_rates)

    seq_diffs = {i:get_diffs(string_seqs[i], hyper_evolved_sequences[i]) for i in range(len(string_seqs))}
    full_filename = os.path.join(
        working_dir, 
        "hyperdata.muts",
    )
    with open(full_filename, 'w') as f:
        json.dump(seq_diffs, f)

    for i, hseq in enumerate(hyper_evolved_sequences):
        sequences[i].seq = Seq.Seq(hseq, alphabet=Alphabet.SingleLetterAlphabet())
    full_filename = os.path.join(
        working_dir, 
        "hyperdata.fasta",
    )
    SeqIO.write(sequences, full_filename, 'fasta')
    return full_filename


def run_proviral(sequences, working_dir, out_dir, platform, paired_end, proviral_fraction):
    """
    Perform hypermutation on a set of sequences and generate an NGS dataset for
    non-hypermutated and hypermutated sequences, then sample these two FASTQ files
    to create a new FASTQ file consisting of a fraction of hypermutate and
    non-hypermutated data.

    Args:
        sequences: path to FASTA sequence file
        working_dir: temporary directory for storing output files 
        out_dir: directory in which to store result
        platform: platform (e.g. roche)
        paired_end: produce paired_end data
        proviral_fraction: fraction of proviral data in final dataset

    Returns:
        True 
    """
    print('Using temporary working directory: {}'.format(working_dir))

    bio_sequences = [s for s in SeqIO.parse(sequences, 'fasta')]
    hypermutated_sequence_file = hypermutate_sequences(bio_sequences, working_dir)

    platf = getattr(plat, platform)

    f = art.simulate(sequences, platf, platf.coverage, paired_end, out_dir)
    fastq_file0, sam_file0 = art.simulate(sequences, platf, platf.coverage, paired_end, out_dir)
    fastq_file1, sam_file1 = art.simulate(hypermutated_sequence_file, platf, platf.coverage, paired_end, out_dir)

    if paired_end:
        fastq_file0a, fastq_file0b = fastq_file0
        fastq_file1a, fastq_file1b = fastq_file1
        open_files = []
        for fq in [fastq_file0a, fastq_file0b, fastq_file1a, fastq_file1b]:
            with open(fq, 'r') as f:
                open_files.append([i for i in parse_fastq(f)])
        fq0a, fq0b, fq1a, fq1b = open_files
    else:
        open_files = []
        for fq in [fastq_file0, fastq_file1]:
            with open(fq, 'r') as f:
                open_files.append([i for i in parse_fastq(f)])
        fq0, fq1 = open_files

    n_reads = min([len(i) for i in open_files])
    n_proviral_reads = int(proviral_fraction * n_reads)
    print('Sampling {0} proviral reads out of {1} total reads (r={2:.3f}).'.format(n_proviral_reads, n_reads, n_proviral_reads/float(n_reads)))

    if paired_end:
    
        mixed_fastq_f, mixed_fastq_r = [], []
        while len(mixed_fastq_f) < n_proviral_reads:
            i_f = random.randint(0, len(fq1a) - 1)
            i_r = random.randint(0, len(fq1b) - 1)
            mixed_fastq_f.append(fq1a.pop(i_f))
            mixed_fastq_r.append(fq1b.pop(i_r))
        while len(mixed_fastq_f) < n_reads:
            i_f = random.randint(0, len(fq0a) - 1)
            i_r = random.randint(0, len(fq0b) - 1)
            mixed_fastq_f.append(fq0a.pop(i_f))
            mixed_fastq_r.append(fq0b.pop(i_r))
       
        mixed_fastq_f = ''.join([j for i in mixed_fastq_f for j in i])
        mixed_fastq_r = ''.join([j for i in mixed_fastq_r for j in i])
    
        full_filename_f = os.path.join(
            out_dir, 
            "mixed_hyperdata1.fastq",
        )
        full_filename_r = os.path.join(
            out_dir, 
            "mixed_hyperdata2.fastq",
        )
        with open(full_filename_f, 'w') as f:
            f.write(mixed_fastq_f) 
        with open(full_filename_r, 'w') as f:
            f.write(mixed_fastq_r) 

    else:
        mixed_fastq = []
        while len(mixed_fastq) < n_proviral_reads:
            i = random.randint(0, len(fq1) - 1)
            mixed_fastq.append(fq1.pop(i))
        while len(mixed_fastq) < n_reads:
            i = random.randint(0, len(fq0) - 1)
            mixed_fastq.append(fq0.pop(i))
   
        mixed_fastq = ''.join([j for i in mixed_fastq for j in i])

        full_filename = os.path.join(
            out_dir, 
            "mixed_hyperdata.fastq",
        )
        with open(full_filename, 'w') as f:
            f.write(mixed_fastq) 

    print 'Output saved in:', out_dir
    return True
