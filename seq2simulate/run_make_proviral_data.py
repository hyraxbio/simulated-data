from Bio import SeqIO, Seq, Alphabet
import os
import random
import json
from glob import glob
import pickle

from hypermutation import hypermutate
import platform as plat
import art
from custom_generators import parse_fastq, parse_sam
import diversity

MUTATIONS = {
    'null': 0,
    'hypermutation': 1,
    'longdel': 2,
    'insertion': 3,
    'frameshift': 4,
    'stopcodon': 5,
    'inversion': 6,
}

# proviral hypermutations per hundred bp
# hypermutation_rate = 3

def get_diffs1(seq0, seq1):
    """
    Simply return a list of 1-based different indices between two strings.
    """
    assert len(seq0) == len(seq1)
    return [i for i in range(len(seq0)) if seq0[i] != seq1[i]]


def run_proviral(sequences_path, working_dir, out_dir, platform, paired_end, proviral_fraction, unclean_working=False, hypermutation_rate=3):
    """
    Perform hypermutation on a set of sequences and generate an NGS dataset for
    non-hypermutated and hypermutated sequences, then sample these two FASTQ files
    to create a new FASTQ file consisting of a fraction of hypermutate and
    non-hypermutated data.

    Args:
        sequences_path: path to FASTA sequence file
        working_dir: temporary directory for storing output files 
        out_dir: directory in which to store result
        platform: platform (e.g. roche)
        paired_end: produce paired_end data
        proviral_fraction: fraction of proviral data in final dataset
        unclean_working: do not delete working directory files upon completion
        hypermutation_rate: rate of hypermutation per 100 bp

    Returns:
        True 


    """
    print('Using temporary working directory: {}'.format(working_dir))

    sequences = [s for s in SeqIO.parse(sequences_path, 'fasta')]

    data_files = {'null': sequences_path}

    sequences_strings = diversity._convert_seqs_to_strs(sequences)
    sequences_strings = diversity._simulate_hypermutation(sequences_strings, hypermutation_rate=hypermutation_rate)
    diversity._update_seq_reads_from_strs(sequences, sequences_strings)
    data_files['hypermutation'] = (_write_to_FASTA(sequences, working_dir, '1_hm_data.fasta'))

    sequences_strings = diversity._convert_seqs_to_strs(sequences)
    sequences_strings = diversity._simulate_deletions(sequences_strings, freq=1)
    diversity._update_seq_reads_from_strs(sequences, sequences_strings)
    data_files['longdel'] = (_write_to_FASTA(sequences, working_dir, '2_del_data.fasta'))

    sequences_strings = diversity._convert_seqs_to_strs(sequences)
    sequences_strings = diversity._simulate_insertions(sequences_strings, freq=1)
    diversity._update_seq_reads_from_strs(sequences, sequences_strings)
    data_files['insertion'] = (_write_to_FASTA(sequences, working_dir, '3_ins_data.fasta'))

    sequences_strings = diversity._convert_seqs_to_strs(sequences)
    sequences_strings = diversity._simulate_frameshifts(sequences_strings, freq=1)
    diversity._update_seq_reads_from_strs(sequences, sequences_strings)
    data_files['frameshift'] = (_write_to_FASTA(sequences, working_dir, '4_fs_data.fasta'))

    sequences_strings = diversity._convert_seqs_to_strs(sequences)
    sequences_strings = diversity._simulate_stop_codons(sequences_strings, freq=1)
    diversity._update_seq_reads_from_strs(sequences, sequences_strings)
    data_files['stopcodon'] = (_write_to_FASTA(sequences, working_dir, '5_sc_data.fasta'))

    sequences_strings = diversity._convert_seqs_to_strs(sequences)
    sequences_strings = diversity._simulate_inversions(sequences_strings, freq=1)
    diversity._update_seq_reads_from_strs(sequences, sequences_strings)
    data_files['inversion'] = (_write_to_FASTA(sequences, working_dir, '6_inv_data.fasta'))


    platf = getattr(plat, platform)

    fastq_files, sam_files = {}, {}
    for mutation_type, data_file in data_files.iteritems():
        fastq_file, sam_file = art.simulate(data_file, platf, platf.coverage, paired_end, working_dir)
        if not isinstance(fastq_file, tuple):
            fastq_file = [fastq_file]
        if not isinstance(sam_file, tuple):
            sam_file = [sam_file]
        fastq_files[mutation_type], sam_files[mutation_type] = [i for i in fastq_file], [i for i in sam_file]

    open_fq_files = {}
    for mutation_type, fq in fastq_files.iteritems():
        fastq_siblings = []
        for _fq in fq:
            with open(_fq, 'r') as f:
                fastq_siblings.append([i for i in parse_fastq(f)])
        open_fq_files[mutation_type] = fastq_siblings
    
    open_sam_files = {}
    for mutation_type, fs in sam_files.iteritems():
        sam_siblings = []
        for _fs in fs:
            with open(_fs, 'r') as f:
                sam_siblings.append([i for i in parse_sam(f)])
        open_sam_files[mutation_type] = sam_siblings


    n_reads = min([len(j) for i in open_fq_files.values() for j in i])
    n_proviral_reads = int(proviral_fraction * n_reads)
    n_reads_frac = int(1/7.0 * n_reads)
    n_mutated_reads = {
        'null': n_reads_frac,
        'hypermutation': n_reads_frac,
        'longdel': n_reads_frac,
        'insertion': n_reads_frac,
        'frameshift': n_reads_frac,
        'stopcodon': n_reads_frac,
        'inversion': n_reads_frac,
    }

    for i,j in n_mutated_reads.iteritems():
        print('Sampling {0} ({1}) reads out of {2} total reads (r={3:.3f}).'.format(j, i, n_reads, j/float(n_reads)))

    # NOTE: 
    # roche: fq sequence is fwd strand, sam seq is fwd strand
    # illumina: fq sequence is rev comp of fwd strand, sam seq is fwd strand
    #           i.e. fq has complement of 3'->5' end of fwd strand, sam has complement of this
    # illumina+paired: 
    #           fq1 sequence is rev comp of fwd strand, sam seq is fwd strand
    #           i.e. fq has complement of 3'->5' end of fwd strand, sam has complement of this
    #           fq2 sequence is fwd strand, sam seq is fwd strand

    fastq_samples = {}
    for mutation_type, fastq_files in open_fq_files.iteritems():
        if len(fastq_files) == 1:
            fastq_samples[mutation_type] = sample_fastq(n_mutated_reads[mutation_type], 
                                                        mutation_code=MUTATIONS[mutation_type], 
                                                        fastq1=fastq_files[0], 
                                                        fastq2=None
                                                       )
        elif len(fastq_files) == 2:
            fastq_samples[mutation_type] = sample_fastq(n_mutated_reads[mutation_type], 
                                                        mutation_code=MUTATIONS[mutation_type], 
                                                        fastq1=fastq_files[0], 
                                                        fastq2=fastq_files[1]
                                                       )
        else:
            raise ValueError('Too many FASTQ files loaded.')

    outfiles = []
    if paired_end:
        mixed_f, mixed_r = [], []
        for fastq_sample in fastq_samples.values(): 
            mixed_f.append(''.join(fastq_sample[0]))
            mixed_r.append(''.join(fastq_sample[1]))
        mixed_f = ''.join(mixed_f)
        mixed_r = ''.join(mixed_r)
        outfiles.append(_write_to_file(mixed_f, out_dir, 'mixed_hyperdata1.fq')) 
        outfiles.append(_write_to_file(mixed_r, out_dir, 'mixed_hyperdata2.fq'))
    else:
        mixed_f = []
        for fastq_sample in fastq_samples.values(): 
            mixed_f.append(''.join(fastq_sample))
        mixed_f = ''.join(mixed_f)
        outfiles.append(_write_to_file(mixed_f, out_dir, 'mixed_hyperdata1.fq')) 

    if unclean_working:
        pass
        # # this is saved primarily for functional testing
        # sam_filename = "sam_reads.pkl"
        # full_filename = os.path.join(
        #     working_dir, 
        #     sam_filename,
        # )
        # with open(full_filename, 'wb') as f:
        #     pickle.dump(sam_reads, f)            
    else:

        for tempfile in glob(working_dir+'/*.fq'):
            os.unlink(tempfile)
        for tempfile in glob(working_dir+'/*.sam'):
            os.unlink(tempfile)
        for tempfile in glob(working_dir+'/*.pkl'):
            os.unlink(tempfile)

    print 'Output saved in:', out_dir
    print('FASTQ mutation codes:')
    for i in range(len(MUTATIONS)):
        for mut, cod in MUTATIONS.iteritems():
            if cod == i:
                print('{}: {}'.format(cod, mut))

    return True

def sample_fastq(n_reads, mutation_code=0, fastq1=None, fastq2=None):

    """
    Returns a random sample of a list form of a FASTQ file (output of
    custom_generators.parse_fastq())

    Args:
        n_reads: number of reads to sample
        mutation_code: see MUTATIONS, this code gets appended to the read ID 
        fastq1/2: open FASTQ file objects (paired-end data has forward and
        reverse files)
    """
    paired_end = False
    if fastq1 is None and fastq2 is None:
        raise TypeError('At least one of fastq1 and fastq2 must not be None.')
    elif fastq1 is None and fastq2 is not None:
        fastq1, fastq2 = fastq2, fastq1
    elif fastq1 is not None and fastq2 is not None:
        paired_end = True

    if paired_end:
        if len(fastq1) != len(fastq2):
            raise ValueError('FASTQ forward and reverse files must have the same lengths.')

    fastq_sample1 = []
    if paired_end:
        fastq_sample2 = []

    while len(fastq_sample1) < n_reads:
        i = random.randint(0, len(fastq1) - 1)
        proviral_read1 = list(fastq1.pop(i))
        read_id1 = proviral_read1[0][1:].strip('\n')
        proviral_read1[0] = proviral_read1[0][:-1] + '_{}\n'.format(mutation_code)
        fastq_sample1.append(''.join(proviral_read1))
        if paired_end:
            proviral_read2 = list(fastq2.pop(i))
            read_id2 = proviral_read2[0][1:].strip('\n')
            proviral_read2[0] = proviral_read2[0][:-1] + '_{}\n'.format(mutation_code)
            fastq_sample2.append(''.join(proviral_read2))

    if paired_end:    
        return fastq_sample1, fastq_sample2
    else:
        return fastq_sample1

def _parse_sam_line(read_id, sam_file, paired_end=False):
    """
    Parse a single SAM file line (or pair of lines for paired-end data).

    Args:
        read_id: id string of FASTQ read
        sam_file: output of custom_generators.parse_sam() 
        paired_end: is this SAM file for paired-end data
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
     

def _write_to_FASTA(sequences, working_dir, filename):
    """
    write sequences to FASTA

    Args:
        sequences: list of BioPython sequence reads
        working_dir: directory for file
        filename: name of file
    
    Returns:
        full path to file
    """
    sequence_file = os.path.join(
        working_dir, 
        filename,
    )
    SeqIO.write(sequences, sequence_file, 'fasta')
    return sequence_file

def _write_to_file(obj, path, filename):
    """
    write object to file

    Args:
        obj: Python object
        working_dir: directory for file
        filename: name of file
    
    Returns:
        full path to file
    """
    full_filename = os.path.join(
        path, 
        filename,
    )
    with open(full_filename, 'w') as f:
        f.write(obj)
    return full_filename
