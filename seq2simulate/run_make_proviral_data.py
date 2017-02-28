from Bio import SeqIO, Seq, Alphabet
import os
import random
import json
from glob import glob
import pickle
import uuid

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


def run_proviral(sequences_path, working_dir, out_dir, platform, paired_end, unclean_working=False, hypermutation_rate=3):
    """

    Perform proviral mutations (including hypermutations, deletions,
    insertions, frameshifts, inserting stop-codons, and inversions) on a set of
    sequences and generate an NGS dataset for viral and proviral sequences, then
    sample these two FASTQ files to create a new FASTQ file consisting of a
    fraction of reads produced from each proviral (and viral) dataset.

    Args:
        sequences_path: path to FASTA sequence file
        working_dir: temporary directory for storing output files 
        out_dir: directory in which to store result
        platform: platform (e.g. roche)
        paired_end: produce paired_end data
        unclean_working: do not delete working directory files upon completion
        hypermutation_rate: rate of hypermutation per 100 bp

    Returns:
        True 


    """
    print('Using temporary working directory: {}'.format(working_dir))

    sequences = [s for s in SeqIO.parse(sequences_path, 'fasta')]

    data_files, diff_files = _make_mutation_data_files(sequences, working_dir, hypermutation_rate=hypermutation_rate)
    data_files['null'] = sequences_path

    fastq_files, sam_files = _make_art_files(data_files, platform, working_dir, paired_end=paired_end)

    open_fq_files, open_sam_files, open_diff_files = _open_all_data_files(sequences, fastq_files, sam_files, diff_files)

    n_reads = min([len(j) for i in open_fq_files.values() for j in i])
    n_reads_frac = int(1.0/len(MUTATIONS) * n_reads)
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
            fastq_files += [None]
        fastq_sample = sample_fastq(n_mutated_reads[mutation_type], 
                                                    fastq1=fastq_files[0], 
                                                    fastq2=fastq_files[1],
                                                   )
        _decorate_fastq_headers(fastq_sample, 
                                mutation_type, 
                                sam_file=open_sam_files[mutation_type], 
                                diff_file=open_diff_files[mutation_type], 
                                paired_end=paired_end)
        fastq_samples[mutation_type] = fastq_sample

    outfiles = []
    output_filestring = str(uuid.uuid4())
    if paired_end:
        mixed_f, mixed_r = [], []
        for fastq_sample in fastq_samples.values(): 
            mixed_f.append(''.join(fastq_sample[0]))
            mixed_r.append(''.join(fastq_sample[1]))
        mixed_f = ''.join(mixed_f)
        mixed_r = ''.join(mixed_r)
        outfiles.append(_write_to_file(mixed_f, out_dir, output_filestring+'_1.fq')) 
        outfiles.append(_write_to_file(mixed_r, out_dir, output_filestring+'_2.fq'))
    else:
        mixed_f = []
        for fastq_sample in fastq_samples.values(): 
            mixed_f.append(''.join(fastq_sample))
        mixed_f = ''.join(mixed_f)
        outfiles.append(_write_to_file(mixed_f, out_dir, output_filestring+'_1.fq')) 

    if not unclean_working:
        for extension in ['fasta', 'data', 'fq', 'sam', 'pkl']:
            for tempfile in glob(working_dir + '/*.' + extension):
                os.unlink(tempfile)

    print 'Output saved in:'
    for outfile in outfiles:
        print(outfile)

    print('FASTQ mutation codes:')
    for i in range(len(MUTATIONS)):
        for mut, cod in MUTATIONS.iteritems():
            if cod == i:
                print('{}: {}'.format(cod, mut))

    return True

def _make_mutation_data_files(sequences, working_dir, hypermutation_rate=3):

    """
    Performs mutations of various kinds on template sequence(s) and saves both
    the mutated sequences and a JSON of the mutation locations to file. Returns
    these files in two dictionaries.
    """
    data_files = {}
    diff_files = {}

    mutation_types  = ['hypermutation',
                       'longdel', 
                       'insertion', 
                       'frameshift', 
                       'stopcodon', 
                       'inversion']

    mutation_funcs = {'hypermutation': [diversity._simulate_hypermutation, {'hypermutation_rate': hypermutation_rate}],
                      'longdel': [diversity._simulate_deletions, {'freq': 1, 'no_frameshifts': True}], 
                      'insertion': [diversity._simulate_insertions, {'freq': 1, 'no_frameshifts': True}], 
                      'frameshift': [diversity._simulate_frameshifts, {'freq': 1}],
                      'stopcodon': [diversity._simulate_stop_codons, {'freq': 1}],
                      'inversion': [diversity._simulate_inversions, {'freq': 1}],
                     }

    for i_mutation, mutation_type in enumerate(mutation_types):
        mutation_func, mutation_kwargs = mutation_funcs[mutation_type]
        sequences_strings, sequence_ids = diversity._convert_seqs_to_strs(sequences)
        sequences_strings, seq_diffs = mutation_func(sequences_strings, **mutation_kwargs)
        seq_diffs = dict(zip(sequence_ids, seq_diffs))
        diversity._update_seq_reads_from_strs(sequences, sequences_strings)
        data_files[mutation_type] = (_write_to_FASTA(sequences, working_dir, '{}_data.fasta'.format(i_mutation)))
        diff_files[mutation_type] = (_write_to_file(seq_diffs, working_dir, '{}_diffs.data'.format(i_mutation)))

    return data_files, diff_files


def _make_art_files(data_files, platform, working_dir, paired_end=False):
    platf = getattr(plat, platform)

    fastq_files, sam_files = {}, {}
    for mutation_type, data_file in data_files.iteritems():
        fastq_file, sam_file = art.simulate(data_file, platf, platf.coverage, paired_end, working_dir)
        if not isinstance(fastq_file, tuple):
            fastq_file = [fastq_file]
        fastq_files[mutation_type], sam_files[mutation_type] = [i for i in fastq_file], sam_file
    return fastq_files, sam_files


def _open_all_data_files(sequences, fastq_files, sam_files, diff_files):
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
    return open_fq_files, open_sam_files, open_diff_files


def _decorate_fastq_headers(fastq_sample, mutation_type, sam_file=None, diff_file=None, paired_end=False):

    """
    Decorate FASTQ headers with a suffix code indicating whether or not this
    read covers a known mutation (see MUTATIONS for mutation codes).
    """
    if paired_end:
        for i_read in range(len(fastq_sample[0])):
            reads = [fastq_sample[0][i_read], fastq_sample[1][i_read]]
            id_suffixes = [None, None]
            if mutation_type != 'null':
                if sam_file is not None and diff_file is not None:
                    read_id = reads[0][0][1:].strip('\n')
                    read_id, read_number = read_id[:-2], read_id[-1]
                    read_seq1, read_seq2 = reads[0][1].strip('\n'), reads[1][1].strip('\n')
                    sam_line_f, sam_line_r = _parse_sam_line(read_id, sam_file, paired_end=paired_end)
                    seq_similarities = [diversity._sim_score(sam_line_f['seq'], read_seq) for read_seq in [read_seq1, read_seq2]]
                    if seq_similarities[1] > seq_similarities[0]:
                        sam_line_f, sam_line_r = sam_line_r, sam_line_f
                    elif seq_similarities[1] == seq_similarities[0]:
                        raise ValueError('Paired-end FASTQ reads could not be distinguished between.')
                    seq_diffs = diff_file[sam_line_f['seq_id']]
                    id_suffixes = [_get_id_suffix(mutation_type, sam_line, seq_diffs) for sam_line in [sam_line_f, sam_line_r]]
            for i in range(len(id_suffixes)):
                if id_suffixes[i] is None:
                    id_suffixes[i] = '_{}\n'.format(MUTATIONS['null'])
                
            reads[0][0] = reads[0][0][:-1] + id_suffixes[0]
            reads[1][0] = reads[1][0][:-1] + id_suffixes[1]

            fastq_sample[0][i_read] = ''.join(reads[0])
            fastq_sample[1][i_read] = ''.join(reads[1])
    else:
        for i_read, read in enumerate(fastq_sample):
            id_suffix = None
            if mutation_type != 'null':
                if sam_file is not None and diff_file is not None:
                    read_id = read[0][1:].strip('\n')

                    sam_line = _parse_sam_line(read_id, sam_file, paired_end=paired_end)
                    seq_diffs = diff_file[sam_line['seq_id']]
                    id_suffix = _get_id_suffix(mutation_type, sam_line, seq_diffs)
            if id_suffix is None:
                id_suffix = '_{}\n'.format(MUTATIONS['null'])
                
            read[0] = read[0][:-1] + id_suffix
            fastq_sample[i_read] = ''.join(read)

def _get_id_suffix(mutation_type, sam_line, seq_diffs):

    mutation_funcs = {'hypermutation': _is_hypermutated,
                      'longdel': _is_longdel,
                      'insertion': _is_insertion,
                      'frameshift': _is_frameshift,
                      'stopcodon': _is_stopcodon,
                      'inversion': _is_inversion,
                     }

    id_suffix = None
    for mutation, is_mutated in mutation_funcs.iteritems():
        if mutation_type == mutation and is_mutated(sam_line, seq_diffs):
            id_suffix = '_{}\n'.format(MUTATIONS[mutation_type])
            break

    return id_suffix 

def _is_hypermutated(sam_line, seq_diffs):
    critical_num_hypermutations = 6
    covered_mutations = [i for i in seq_diffs if i >= sam_line['read_start'] and i <= sam_line['read_end']]
    if len(covered_mutations) >= critical_num_hypermutations:
        return True
    return False

def _is_longdel(sam_line, seq_diffs):
    if sam_line['read_start'] <= seq_diffs[0] and sam_line['read_end'] > seq_diffs[1]:
        return True
    return False

def _is_insertion(sam_line, seq_diffs):
    if seq_diffs[0] <= sam_line['read_start'] and seq_diffs[1] > sam_line['read_start'] \
        or sam_line['read_start'] <= seq_diffs[0] <= sam_line['read_end']:
        return True
    return False

def _is_frameshift(sam_line, seq_diffs):
    if sam_line['read_start'] < seq_diffs[0] and sam_line['read_end'] > seq_diffs[0] \
        or sam_line['read_start'] >= seq_diffs[0] and sam_line['read_end'] > seq_diffs[0]:
        return True
    return False

def _is_stopcodon(sam_line, seq_diffs):
    if sam_line['read_start'] <= seq_diffs[0] and sam_line['read_end'] >= seq_diffs[0]:
        return True
    return False

def _is_inversion(sam_line, seq_diffs):
    if seq_diffs[0] <= sam_line['read_start'] and seq_diffs[1] > sam_line['read_start'] \
        or sam_line['read_start'] <= seq_diffs[0] <= sam_line['read_end']:
        return True
    return False

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
        return [fastq_sample1, fastq_sample2]
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
    seq_ind = 8
    if paired_end:
        seq_id_1 = sam_file[read_id][0][1] 
        seq_id_2 = sam_file[read_id][1][1] 
        if seq_id_1 != seq_id_2:
            raise ValueError('read_id does not reference paired reads from same sequence: {} {}'.format(seq_id_1, seq_id_2))
        # get the FLAG bit for directionality 
        if (int(sam_file[read_id][0][0]) & 0x10) and not (int(sam_file[read_id][1][0]) & 0x10):
            sam_file[read_id] = sam_file[read_id][::-1]
        elif int(sam_file[read_id][0][0]) == int(sam_file[read_id][1][0]):
            raise ValueError('Both paired-end FASTQ reads are in the same direction.')
        
        read_start_f = int(sam_file[read_id][0][2])
        read_start_r = int(sam_file[read_id][1][2])
        read_end_f = read_start_f+len(sam_file[read_id][0][seq_ind])-1
        read_end_r = read_start_r+len(sam_file[read_id][1][seq_ind])-1

        if read_start_f >= read_end_r:
            raise ValueError('Malformed FASTQ. Paired-end read directions incorrectly specified.')

        result = [
                  {'seq_id':seq_id_1, 'read_start': read_start_f, 'read_end':read_end_f, 'seq':sam_file[read_id][0][seq_ind]},
                  {'seq_id':seq_id_1, 'read_start': read_start_r, 'read_end':read_end_r, 'seq':sam_file[read_id][1][seq_ind]},
                 ] 
    else:
        seq_id = sam_file[read_id][0][1] 
        read_start = int(sam_file[read_id][0][2])
        read_end = read_start+len(sam_file[read_id][0][seq_ind])-1
        result = {'seq_id':seq_id, 'read_start': read_start, 'read_end':read_end, 'seq':sam_file[read_id][0][seq_ind]}
    return result
     

def _write_to_FASTA(sequences, working_dir, filename):
    """
    Write sequences to FASTA.

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
    Write object to file. If not str, pickle.

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

 
