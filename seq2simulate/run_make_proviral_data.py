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
    diff_files = {}

    sequences_strings, sequence_ids = diversity._convert_seqs_to_strs(sequences)
    sequences_strings, seq_diffs = diversity._simulate_hypermutation(sequences_strings, hypermutation_rate=hypermutation_rate)
    seq_diffs = dict(zip(sequence_ids, seq_diffs))
    diversity._update_seq_reads_from_strs(sequences, sequences_strings)
    data_files['hypermutation'] = (_write_to_FASTA(sequences, working_dir, '1_hm_data.fasta'))
    diff_files['hypermutation'] = (_write_to_file(seq_diffs, working_dir, '1_hm_diffs.data'))

    sequences_strings, sequence_ids = diversity._convert_seqs_to_strs(sequences)
    sequences_strings, seq_diffs = diversity._simulate_deletions(sequences_strings, freq=1)
    seq_diffs = dict(zip(sequence_ids, seq_diffs))
    diversity._update_seq_reads_from_strs(sequences, sequences_strings)
    data_files['longdel'] = (_write_to_FASTA(sequences, working_dir, '2_del_data.fasta'))
    diff_files['longdel'] = (_write_to_file(seq_diffs, working_dir, '2_del_diffs.data'))

    sequences_strings, sequence_ids = diversity._convert_seqs_to_strs(sequences)
    sequences_strings, seq_diffs = diversity._simulate_insertions(sequences_strings, freq=1)
    seq_diffs = dict(zip(sequence_ids, seq_diffs))
    diversity._update_seq_reads_from_strs(sequences, sequences_strings)
    data_files['insertion'] = (_write_to_FASTA(sequences, working_dir, '3_ins_data.fasta'))
    diff_files['insertion'] = (_write_to_file(seq_diffs, working_dir, '3_ins_diffs.data'))

    sequences_strings, sequence_ids = diversity._convert_seqs_to_strs(sequences)
    sequences_strings, seq_diffs = diversity._simulate_frameshifts(sequences_strings, freq=1)
    seq_diffs = dict(zip(sequence_ids, seq_diffs))
    diversity._update_seq_reads_from_strs(sequences, sequences_strings)
    data_files['frameshift'] = (_write_to_FASTA(sequences, working_dir, '4_fs_data.fasta'))
    diff_files['frameshift'] = (_write_to_file(seq_diffs, working_dir, '4_fs_diffs.data'))

    sequences_strings, sequence_ids = diversity._convert_seqs_to_strs(sequences)
    sequences_strings, seq_diffs = diversity._simulate_stop_codons(sequences_strings, freq=1)
    seq_diffs = dict(zip(sequence_ids, seq_diffs))
    diversity._update_seq_reads_from_strs(sequences, sequences_strings)
    data_files['stopcodon'] = (_write_to_FASTA(sequences, working_dir, '5_sc_data.fasta'))
    diff_files['stopcodon'] = (_write_to_file(seq_diffs, working_dir, '5_sc_diffs.data'))

    sequences_strings, sequence_ids = diversity._convert_seqs_to_strs(sequences)
    sequences_strings, seq_diffs = diversity._simulate_inversions(sequences_strings, freq=1)
    seq_diffs = dict(zip(sequence_ids, seq_diffs))
    diversity._update_seq_reads_from_strs(sequences, sequences_strings)
    data_files['inversion'] = (_write_to_FASTA(sequences, working_dir, '6_inv_data.fasta'))
    diff_files['inversion'] = (_write_to_file(seq_diffs, working_dir, '5_inv_diffs.data'))


    platf = getattr(plat, platform)

    fastq_files, sam_files = {}, {}
    for mutation_type, data_file in data_files.iteritems():
        fastq_file, sam_file = art.simulate(data_file, platf, platf.coverage, paired_end, working_dir)
        if not isinstance(fastq_file, tuple):
            fastq_file = [fastq_file]
        fastq_files[mutation_type], sam_files[mutation_type] = [i for i in fastq_file], sam_file

    open_fq_files = {}
    for mutation_type, fq in fastq_files.iteritems():
        fastq_siblings = []
        for _fq in fq:
            with open(_fq, 'r') as f:
                fastq_siblings.append([i for i in parse_fastq(f)])
        open_fq_files[mutation_type] = fastq_siblings
    
    open_sam_files = {}
    for mutation_type, fs in sam_files.iteritems():
        with open(fs, 'r') as f:
            sam_reads = {}
            for sam_read in parse_sam(f):
                read_id = sam_read[0].strip('\n') 
                if read_id in sam_reads:
                    sam_reads[read_id].append(sam_read[1:])
                else:
                    sam_reads[read_id] = [sam_read[1:]]
            open_sam_files[mutation_type] = sam_reads

    open_diff_files = {'null': {i.id: [] for i in sequences}}
    for mutation_type, fd in diff_files.iteritems():
        with open(fd, 'r') as f:
            open_diff_files[mutation_type] = pickle.load(f)

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
            fastq_sample = sample_fastq(n_mutated_reads[mutation_type], 
                                                        fastq1=fastq_files[0], 
                                                        fastq2=None
                                                       )
            _decorate_fastq_headers(fastq_sample, 
                                    mutation_type, 
                                    sam_file=open_sam_files[mutation_type], 
                                    diff_file=open_diff_files[mutation_type], 
                                    paired_end=False)
            fastq_samples[mutation_type] = fastq_sample
        elif len(fastq_files) == 2:
            fastq_sample1, fastq_sample2 = sample_fastq(n_mutated_reads[mutation_type], 
                                                        fastq1=fastq_files[0], 
                                                        fastq2=fastq_files[1]
                                                       )
            _decorate_fastq_headers(fastq_sample1, 
                                    mutation_type, 
                                    sam_file=open_sam_files[mutation_type], 
                                    diff_file=open_diff_files[mutation_type], 
                                    paired_end=True)
            _decorate_fastq_headers(fastq_sample2, 
                                    mutation_type, 
                                    sam_file=open_sam_files[mutation_type], 
                                    diff_file=open_diff_files[mutation_type], 
                                    paired_end=True)
            fastq_samples[mutation_type] = [fastq_sample1, fastq_sample2]
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

def _decorate_fastq_headers(fastq_sample, mutation_type, sam_file=None, diff_file=None, paired_end=False):

    critical_num_hypermutations = 6

    for i_read, read in enumerate(fastq_sample):
        id_suffix = '_{}\n'.format(MUTATIONS['null'])
        if sam_file is not None and diff_file is not None:
            read_id = read[0][1:].strip('\n')
            if paired_end:
                read_id, read_direction = read_id[:-2], read_id[-1]
                sam_line_f, sam_line_r = _parse_sam_line(read_id, sam_file, paired_end=paired_end)
                if read_direction == '1':
                    sam_line = sam_line_f
                elif read_direction == '2':
                    sam_line = sam_line_r
                else:
                    raise ValueError('Paired-end FASTQ read headers must have directionality expressed as .../1 or .../2')
            else:
                sam_line = _parse_sam_line(read_id, sam_file, paired_end=paired_end)
            seq_diffs = diff_file[sam_line['seq_id']]

            if mutation_type == 'hypermutation':
                covered_mutations = [i for i in seq_diffs if i >= sam_line['read_start'] and i <= sam_line['read_end']]
                if len(covered_mutations) >= critical_num_hypermutations:
                    id_suffix = '_{}\n'.format(MUTATIONS[mutation_type])
            elif mutation_type == 'longdel':
                if sam_line['read_start'] <= seq_diffs[0] and sam_line['read_end'] > seq_diffs[1]:
                    id_suffix = '_{}\n'.format(MUTATIONS[mutation_type])
            elif mutation_type == 'insertion':
                if seq_diffs[0] <= sam_line['read_start'] <= seq_diffs[1] \
                    or seq_diffs[0] <= sam_line['read_end'] <= seq_diffs[1] \
                    or sam_line['read_start'] <= seq_diffs[0] and sam_line['read_end'] >= seq_diffs[1] \
                    or sam_line['read_start'] >= seq_diffs[0] and sam_line['read_end'] <= seq_diffs[1]:
                    id_suffix = '_{}\n'.format(MUTATIONS[mutation_type])
            elif mutation_type == 'frameshift':
                if sam_line['read_start'] < seq_diffs[0] and sam_line['read_end'] > seq_diffs[0] \
                    or sam_line['read_start'] >= seq_diffs[0] and sam_line['read_end'] > seq_diffs[0]:
                    id_suffix = '_{}\n'.format(MUTATIONS[mutation_type])
            elif mutation_type == 'stopcodon':
                if sam_line['read_start'] <= seq_diffs[0] and sam_line['read_end'] >= seq_diffs[0]:
                    id_suffix = '_{}\n'.format(MUTATIONS[mutation_type])
            elif mutation_type == 'inversion':
                if seq_diffs[0] <= sam_line['read_start'] <= seq_diffs[1] \
                    or seq_diffs[0] <= sam_line['read_end'] <= seq_diffs[1] \
                    or sam_line['read_start'] <= seq_diffs[0] and sam_line['read_end'] >= seq_diffs[1] \
                    or sam_line['read_start'] >= seq_diffs[0] and sam_line['read_end'] <= seq_diffs[1]:
                    id_suffix = '_{}\n'.format(MUTATIONS[mutation_type])
                

            """continue here by getting the read coverage from sam_line and
            checking if it covers the diffs in diff_file"""
        read[0] = read[0][:-1] + id_suffix

        fastq_sample[i_read] = ''.join(read)

def sample_fastq(n_reads, fastq1=None, fastq2=None):

    """
    Returns a random sample of a list form of a FASTQ file (output of
    custom_generators.parse_fastq())

    Args:
        n_reads: number of reads to sample
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
        fastq_sample1.append(proviral_read1)
        if paired_end:
            proviral_read2 = list(fastq2.pop(i))
            fastq_sample2.append(proviral_read2)

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
        seq_id_1 = sam_file[read_id][0][1] 
        seq_id_2 = sam_file[read_id][1][1] 
        if seq_id_1 != seq_id_2:
            raise ValueError('read_id does not reference paired reads from same sequence: {} {}'.format(seq_id_1, seq_id_2))
        if int(sam_file[read_id][0][7]) < 0 and int(sam_file[read_id][1][7]) > 0:
            sam_file[read_id] = sam_file[read_id][::-1]
        elif int(sam_file[read_id][0][7]) > 0 and int(sam_file[read_id][1][7]) > 0 or int(sam_file[read_id][0][7]) < 0 and int(sam_file[read_id][1][7]) < 0:
            raise ValueError('Both paired-end FASTQ reads are in the same direction.')
        
        read_start_f = int(sam_file[read_id][0][2])
        read_start_r = int(sam_file[read_id][1][2])
        seq_ind = 9
        read_end_f = read_start_f+len(sam_file[read_id][0][seq_ind])-1
        read_end_r = read_start_r+len(sam_file[read_id][1][seq_ind])-1

        if read_start_f >= read_end_r:
            raise ValueError('Malformed FASTQ. Paired-end read directions incorrectly specified.')

        result = [
                  {'seq_id':seq_id_1, 'read_start': read_start_f, 'read_end':read_end_f},
                  {'seq_id':seq_id_1, 'read_start': read_start_r, 'read_end':read_end_r},
                 ] 
    else:
        seq_id = sam_file[read_id][0][1] 
        read_start = int(sam_file[read_id][0][2])
        seq_ind = 8
        read_end = read_start+len(sam_file[read_id][0][seq_ind])-1
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
        if isinstance(obj, str):
            f.write(obj)
        else:
            pickle.dump(obj, f)
    return full_filename

 
