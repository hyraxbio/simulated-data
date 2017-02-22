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

def hypermutate_sequences(sequences, hypermutation_rate=3):
    """
    Hypermutates a set of sequence strings.

    Args:
        sequences: list of DNA strings

    Returns:
        list of sequences
        dict of sequence differences
    """
    string_seqs = [str(s.seq) for s in sequences]
    string_seqs_ids = [str(s.id) for s in sequences]
    if len(set(string_seqs_ids)) != len(string_seqs_ids):
        raise ValueError('FASTA sequences must have unique IDs')
    hyper_evolved_sequences = diversity._simulate_hypermutation(string_seqs, hypermutation_rate=hypermutation_rate)

    seq_diffs = {string_seqs_ids[i]: get_diffs1(string_seqs[i], hyper_evolved_sequences[i]) for i in range(len(string_seqs))}

    for i, hseq in enumerate(hyper_evolved_sequences):
        sequences[i].seq = Seq.Seq(hseq, alphabet=Alphabet.SingleLetterAlphabet())

    return sequences, seq_diffs


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
    # hypermutated_sequences, hypermutated_diffs = hypermutate_sequences(bio_sequences, hypermutation_rate=hypermutation_rate)

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

    
   
    # # write sequences to FASTA
    # hypermutated_sequence_file = os.path.join(
    #     working_dir, 
    #     "hyperdata.fasta",
    # )
    # SeqIO.write(sequences, hypermutated_sequence_file, 'fasta')

    # # write JSON of hypermutations
    # hypermutated_diffs_filename = os.path.join(
    #     out_dir, 
    #     "hyperdata.muts",
    # )
    # with open(hypermutated_diffs_filename, 'w') as f:
    #     json.dump(hypermutated_diffs, f)


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

    print('Continue from here...')
    assert False

    # if paired_end:
    #     """
    #     viral    fwd: fq0a, rev: fq0b
    #     proviral fwd: fq1a, rev: fq1b
    #     """
    #     mixed_fastq_f, mixed_fastq_r = [], []
    #     sam_reads = {}
    #     while len(mixed_fastq_f) < n_proviral_reads:
    #         i = random.randint(0, len(fq1a) - 1)

    #         proviral_read_f = list(fq1a.pop(i))
    #         proviral_read_r = list(fq1b.pop(i))
    #         read_id_f = proviral_read_f[0][1:].strip('\n')
    #         read_id_r = proviral_read_r[0][1:].strip('\n')
    #         read_id = read_id_f[:-2]
    #         sam_read = _parse_sam_line(read_id, fs1, paired_end=paired_end)
    #         sam_reads[read_id] = sam_read
    #         #n_hypermutations_f, n_hypermutations_r = _get_n_hypermutations(sam_read, hypermutated_diffs, paired_end=paired_end)
    #         #proviral_read_f[0] = proviral_read_f[0][:-1] + '_{}\n'.format(n_hypermutations_f)
    #         #proviral_read_r[0] = proviral_read_r[0][:-1] + '_{}\n'.format(n_hypermutations_r)
    #         proviral_read_f[0] = proviral_read_f[0][:-1] + '_{}\n'.format(MUTATIONS['hypermutation'])
    #         proviral_read_r[0] = proviral_read_r[0][:-1] + '_{}\n'.format(MUTATIONS['hypermutation'])
    #         mixed_fastq_f.append(proviral_read_f)
    #         mixed_fastq_r.append(proviral_read_r)

    #     while len(mixed_fastq_f) < n_reads:
    #         i = random.randint(0, len(fq0a) - 1)

    #         proviral_read_f = list(fq0a.pop(i))
    #         proviral_read_r = list(fq0b.pop(i))
    #         proviral_read_f[0] = proviral_read_f[0][:-1] + '_{}\n'.format(MUTATIONS['null'])
    #         proviral_read_r[0] = proviral_read_r[0][:-1] + '_{}\n'.format(MUTATIONS['null'])
    #         mixed_fastq_f.append(proviral_read_f)
    #         mixed_fastq_r.append(proviral_read_r)
    #    
    #     mixed_fastq_f = ''.join([j for i in mixed_fastq_f for j in i])
    #     mixed_fastq_r = ''.join([j for i in mixed_fastq_r for j in i])
    # 
    #     full_filename_f = os.path.join(
    #         out_dir, 
    #         "mixed_hyperdata1.fq",
    #     )
    #     full_filename_r = os.path.join(
    #         out_dir, 
    #         "mixed_hyperdata2.fq",
    #     )
    #     with open(full_filename_f, 'w') as f:
    #         f.write(mixed_fastq_f) 
    #     with open(full_filename_r, 'w') as f:
    #         f.write(mixed_fastq_r) 

    # else:
    #     # mixed_fastq = []
    #     # #sam_reads = {}
    #     # while len(mixed_fastq) < n_proviral_reads:
    #     #     i = random.randint(0, len(fq1) - 1)
    #     #     proviral_read = list(fq1.pop(i))
    #     #     read_id = proviral_read[0][1:].strip('\n')
    #     #     #sam_read = _parse_sam_line(read_id, fs1, paired_end=paired_end)
    #     #     #sam_reads[read_id] = sam_read
    #     #     # n_hypermutations = _get_n_hypermutations(sam_read, hypermutated_diffs, paired_end=paired_end)
    #     #     # proviral_read[0] = proviral_read[0][:-1] + '_{}\n'.format(n_hypermutations)
    #     #     proviral_read[0] = proviral_read[0][:-1] + '_{}\n'.format(MUTATIONS['hypermutation'])
    #     #     mixed_fastq.append(proviral_read)
    #     # while len(mixed_fastq) < n_reads:
    #     #     i = random.randint(0, len(fq0) - 1)
    #     #     viral_read = list(fq0.pop(i))
    #     #     viral_read[0] = viral_read[0][:-1] + '_{}\n'.format(MUTATIONS['null'])
    #     #     mixed_fastq.append(viral_read)
   
    #     # mixed_fastq = ''.join([j for i in mixed_fastq for j in i])
    #          
    #     short_filename = "mixed_hyperdata.fq"
    #     full_filename = os.path.join(
    #         out_dir, 
    #         short_filename,
    #     )
    #     with open(full_filename, 'w') as f:
    #         f.write(mixed_fastq) 

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
    return True

def sample_fastq(n_reads, mutation_code=0, fastq1=None, fastq2=None):

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

    #sam_reads = {}
    while len(fastq_sample1) < n_reads:
        i = random.randint(0, len(fastq1) - 1)
        proviral_read1 = list(fastq1.pop(i))
        read_id1 = proviral_read1[0][1:].strip('\n')
        proviral_read1[0] = proviral_read1[0][:-1] + '_{}\n'.format(mutation_code)
        fastq_sample1.append(proviral_read1)
        if paired_end:
            proviral_read2 = list(fastq2.pop(i))
            read_id2 = proviral_read2[0][1:].strip('\n')
            proviral_read2[0] = proviral_read2[0][:-1] + '_{}\n'.format(mutation_code)
            fastq_sample2.append(proviral_read2)

        #sam_read = _parse_sam_line(read_id, fs1, paired_end=paired_end)
        #sam_reads[read_id] = sam_read
        # n_hypermutations = _get_n_hypermutations(sam_read, hypermutated_diffs, paired_end=paired_end)
        # proviral_read[0] = proviral_read[0][:-1] + '_{}\n'.format(n_hypermutations)

    if paired_end:    
        return fastq_sample1, fastq_sample2
    else:
        return fastq_sample1

def _parse_sam_line(read_id, sam_file, paired_end=False):
    """
    Parse a single SAM file line (or pair of lines for paired-end data).
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
     

def _get_n_hypermutations(sam_read, hypermutations, paired_end=False):
    """
    For a given read, return number of hypermutations.

    Args:
        sam_read: output of _parse_sam_line 
        hypermutations: dictionary with sequence ids as keys and lists of hypermutations indices as values
        paired_end: is this paired-end data 
    """
    #sam_read = _parse_sam_line(read_id, sam_file, paired_end=paired_end)
    if paired_end:
        result_f = len([i for i in hypermutations[sam_read['seq_id']] if sam_read['read_start_f'] <= i <= sam_read['read_end_f']])
        result_r = len([i for i in hypermutations[sam_read['seq_id']] if sam_read['read_start_r'] <= i <= sam_read['read_end_r']])
        result = [result_f, result_r]
    else:
        result = len([i for i in hypermutations[sam_read['seq_id']] if sam_read['read_start'] <= i <= sam_read['read_end']])
    return result
    

def _write_to_FASTA(sequences, working_dir, filename):
    # write sequences to FASTA
    sequence_file = os.path.join(
        working_dir, 
        filename,
    )
    SeqIO.write(sequences, sequence_file, 'fasta')
    return sequence_file
