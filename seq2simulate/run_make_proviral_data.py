from Bio import SeqIO, Seq, Alphabet
import os
import random
import json

from hypermutation import hypermutate
import platform as plat
import art

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

def parse_fastq(filename):
    """
    Parse a FASTQ file.
    
    Args:
        filename: path to FASTQ file.
    
    Returns:
        A reshaped list representation of the FASTQ file. Each read is a list of four lines.
    """
    try:
        with open(filename, 'r') as f:
            fq = f.read()
        fq = fq.split('\n')
        fq = ['\n'.join(fq[i*4:i*4+4]) for i in range(len(fq)/4)]
        return fq
    except IOError:
        print('File not found.')


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

    fastq_file0, sam_file0 = art.simulate(sequences, platf, platf.coverage, paired_end, out_dir)
    fastq_file1, sam_file1 = art.simulate(hypermutated_sequence_file, platf, platf.coverage, paired_end, out_dir)

    fq0 = parse_fastq(fastq_file0) 
    fq1 = parse_fastq(fastq_file1) 

    n_reads = min(len(fq0), len(fq1))
    n_proviral_reads = int(proviral_fraction * n_reads)

    print('Sampling {0} proviral reads out of {1} total reads (r={2:.3f}).'.format(n_proviral_reads, n_reads, n_proviral_reads/float(n_reads)))

    mixed_fastq = []
    while len(mixed_fastq) < n_proviral_reads:
        i = random.randint(0, len(fq1) - 1)
        mixed_fastq.append(fq1.pop(i))
    while len(mixed_fastq) < n_reads:
        i = random.randint(0, len(fq0) - 1)
        mixed_fastq.append(fq0.pop(i))
   
    mixed_fastq = '\n'.join(mixed_fastq)

    full_filename = os.path.join(
        out_dir, 
        "mixed_hyperdata.fastq",
    )
    with open(full_filename, 'w') as f:
        f.write(mixed_fastq) 

    print 'Output saved in:', out_dir
    return True
