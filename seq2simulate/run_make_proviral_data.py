from Bio import SeqIO, Seq, Alphabet
import os
import random
import json
from glob import glob

from hypermutation import hypermutate
import platform as plat
import art
from custom_generators import parse_fastq, parse_sam

# proviral hypermutations per hundred bp
hypermutation_rate = 3

def get_diffs1(seq0, seq1):
    """
    Simply return a list of 1-based different indices between two strings.
    """
    assert len(seq0) == len(seq1)
    return [i for i in range(len(seq0)) if seq0[i] != seq1[i]]

def hypermutate_sequences(sequences):
    """
    Hypermutates a set of sequence strings.

    Args:
        sequences: list of DNA strings

    Returns:
        FASTA filename of hypermutated sequences
    """
    print('\n-------------------------------------------------------')
    print('Hypermutating evolved sequences (rate = {} per 100 bp).'.format(hypermutation_rate))
    print('-------------------------------------------------------\n')
    string_seqs = [str(s.seq) for s in sequences]
    string_seqs_ids = [str(s.id) for s in sequences]
    if len(set(string_seqs_ids)) != len(string_seqs_ids):
        raise ValueError('FASTA sequences must have unique IDs')
    hypermutation_rates = [int(hypermutation_rate * len(s) // 100) for s in string_seqs]
    hyper_evolved_sequences = hypermutate.mutate_sequences(string_seqs, hypermutation_rates)

    seq_diffs = {string_seqs_ids[i]: get_diffs1(string_seqs[i], hyper_evolved_sequences[i]) for i in range(len(string_seqs))}

    for i, hseq in enumerate(hyper_evolved_sequences):
        sequences[i].seq = Seq.Seq(hseq, alphabet=Alphabet.SingleLetterAlphabet())

    return sequences, seq_diffs


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
    hypermutated_sequences, hypermutated_diffs = hypermutate_sequences(bio_sequences)

    # write sequences to FASTA
    hypermutated_sequence_file = os.path.join(
        working_dir, 
        "hyperdata.fasta",
    )
    SeqIO.write(hypermutated_sequences, hypermutated_sequence_file, 'fasta')

    # write JSON of hypermutations
    hypermutated_diffs_filename = os.path.join(
        out_dir, 
        "hyperdata.muts",
    )
    with open(hypermutated_diffs_filename, 'w') as f:
        json.dump(hypermutated_diffs, f)

    platf = getattr(plat, platform)

    fastq_file0, sam_file0 = art.simulate(sequences, platf, platf.coverage, paired_end, working_dir)
    fastq_file1, sam_file1 = art.simulate(hypermutated_sequence_file, platf, platf.coverage, paired_end, working_dir)


    if paired_end:
        fastq_file0a, fastq_file0b = fastq_file0
        fastq_file1a, fastq_file1b = fastq_file1
        open_fq_files = []
        for fq in [fastq_file0a, fastq_file0b, fastq_file1a, fastq_file1b]:
            with open(fq, 'r') as f:
                open_fq_files.append([i for i in parse_fastq(f)])
        fq0a, fq0b, fq1a, fq1b = open_fq_files
        open_sam_files = []
        for fs in [sam_file0, sam_file1]:
            open_sam_file = {}
            with open(fs, 'r') as f:
                for i in parse_sam(f):
                    if i[0] in open_sam_file:
                        open_sam_file[i[0]].append(i[1:])
                    else:
                        open_sam_file[i[0]] = [i[1:]]
            open_sam_files.append(open_sam_file)
        fs0, fs1 = open_sam_files
    else:
        open_fq_files = []
        for fq in [fastq_file0, fastq_file1]:
            with open(fq, 'r') as f:
                open_fq_files.append([i for i in parse_fastq(f)])
        fq0, fq1 = open_fq_files
        open_sam_files = []
        for fs in [sam_file0, sam_file1]:
            with open(fs, 'r') as f:
                open_sam_files.append({i[0]:i[1:] for i in parse_sam(f)})
        fs0, fs1 = open_sam_files

    n_reads = min([len(i) for i in open_fq_files])
    n_proviral_reads = int(proviral_fraction * n_reads)
    print('Sampling {0} proviral reads out of {1} total reads (r={2:.3f}).'.format(n_proviral_reads, n_reads, n_proviral_reads/float(n_reads)))

    # NOTE: 
    # roche: fq sequence is fwd strand, sam seq is fwd strand
    # illumina: fq sequence is rev comp of fwd strand, sam seq is fwd strand
    #           i.e. fq has complement of 3'->5' end of fwd strand, sam has complement of this
    # illumina+paired: 
    #           fq1 sequence is rev comp of fwd strand, sam seq is fwd strand
    #           i.e. fq has complement of 3'->5' end of fwd strand, sam has complement of this
    #           fq2 sequence is fwd strand, sam seq is fwd strand
    if paired_end:
        """
        viral    fwd: fq0a, rev: fq0b
        proviral fwd: fq1a, rev: fq1b
        """
        mixed_fastq_f, mixed_fastq_r = [], []
        while len(mixed_fastq_f) < n_proviral_reads:
            i = random.randint(0, len(fq1a) - 1)

            proviral_read_f = list(fq1a.pop(i))
            proviral_read_r = list(fq1b.pop(i))
            read_id_f = proviral_read_f[0][1:].strip('\n')
            read_id_r = proviral_read_r[0][1:].strip('\n')
            read_id = read_id_f[:-2]
            n_hypermutations_f, n_hypermutations_r = _get_n_hypermutations(hypermutated_diffs, read_id, fs1, paired_end=paired_end)
            proviral_read_f[0] = proviral_read_f[0][:-1] + '\t{}\n'.format(n_hypermutations_f)
            proviral_read_r[0] = proviral_read_r[0][:-1] + '\t{}\n'.format(n_hypermutations_r)
            mixed_fastq_f.append(proviral_read_f)
            mixed_fastq_r.append(proviral_read_r)

        while len(mixed_fastq_f) < n_reads:
            i = random.randint(0, len(fq0a) - 1)

            proviral_read_f = list(fq0a.pop(i))
            proviral_read_r = list(fq0b.pop(i))
            proviral_read_f[0] = proviral_read_f[0][:-1] + '\t{}\n'.format(0)
            proviral_read_r[0] = proviral_read_r[0][:-1] + '\t{}\n'.format(0)
            mixed_fastq_f.append(proviral_read_f)
            mixed_fastq_r.append(proviral_read_r)
       
        mixed_fastq_f = ''.join([j for i in mixed_fastq_f for j in i])
        mixed_fastq_r = ''.join([j for i in mixed_fastq_r for j in i])
    
        full_filename_f = os.path.join(
            out_dir, 
            "mixed_hyperdata1.fq",
        )
        full_filename_r = os.path.join(
            out_dir, 
            "mixed_hyperdata2.fq",
        )
        with open(full_filename_f, 'w') as f:
            f.write(mixed_fastq_f) 
        with open(full_filename_r, 'w') as f:
            f.write(mixed_fastq_r) 

    else:
        mixed_fastq = []
        while len(mixed_fastq) < n_proviral_reads:
            i = random.randint(0, len(fq1) - 1)
            proviral_read = list(fq1.pop(i))
            read_id = proviral_read[0][1:].strip('\n')
            n_hypermutations = _get_n_hypermutations(hypermutated_diffs, read_id, fs1, paired_end=paired_end)
            proviral_read[0] = proviral_read[0][:-1] + '\t{}\n'.format(n_hypermutations)
            mixed_fastq.append(proviral_read)
        while len(mixed_fastq) < n_reads:
            i = random.randint(0, len(fq0) - 1)
            viral_read = list(fq0.pop(i))
            viral_read[0] = viral_read[0][:-1] + '\t{}\n'.format(0)
            mixed_fastq.append(viral_read)
   
        mixed_fastq = ''.join([j for i in mixed_fastq for j in i])

        short_filename = "mixed_hyperdata.fq"
        full_filename = os.path.join(
            out_dir, 
            short_filename,
        )
        with open(full_filename, 'w') as f:
            f.write(mixed_fastq) 

    for tempfile in glob(working_dir+'/*.fq'):
        os.unlink(tempfile)
    for tempfile in glob(working_dir+'/*.sam'):
        os.unlink(tempfile)

    print 'Output saved in:', out_dir
    return True

def _parse_sam_line(read_id, sam_file, paired_end=False):
    """
    """
    if paired_end:
        seq_id = sam_file[read_id][0][1] 
        read_start_f = int(sam_file[read_id][0][2])
        read_start_r = int(sam_file[read_id][1][2])
        seq_ind = 9
        read_end_f = read_start_f+len(sam_file[read_id][0][seq_ind])-1
        read_end_r = read_start_r+len(sam_file[read_id][1][seq_ind])-1
        result = {'seq_id':seq_id, 
                  'read_start_f': read_start_f, 'read_end_f':read_end_f,
                  'read_start_r': read_start_r, 'read_end_r':read_end_r,
                 }
    else:
        seq_id = sam_file[read_id][1] 
        read_start = int(sam_file[read_id][2])
        seq_ind = 8
        read_end = read_start+len(sam_file[read_id][seq_ind])-1
        result = {'seq_id':seq_id, 'read_start': read_start, 'read_end':read_end}
    return result
     

def _get_n_hypermutations(hypermutations, read_id, sam_file, paired_end=False):
    """
    For a given read, return number of hypermutations.
    """
    sam_read = _parse_sam_line(read_id, sam_file, paired_end=paired_end)
    if paired_end:
        result_f = len([i for i in hypermutations[sam_read['seq_id']] if sam_read['read_start_f'] <= i <= sam_read['read_end_f']])
        result_r = len([i for i in hypermutations[sam_read['seq_id']] if sam_read['read_start_r'] <= i <= sam_read['read_end_r']])
        result = [result_f, result_r]
    else:
        result = len([i for i in hypermutations[sam_read['seq_id']] if sam_read['read_start'] <= i <= sam_read['read_end']])
    return result
    


