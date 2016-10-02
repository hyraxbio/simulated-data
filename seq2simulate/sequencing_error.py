import copy
import os
import random
import shutil
import subprocess
import uuid

import Bio
import pysam

import art

import custom_generators
import platform as plat
import seq2simulate

# values for ion torrent and miseq taken from
# http://aem.asm.org/content/early/2014/09/22/AEM.02206-14.full.pdf+html

simulated_fq_name = "_art.fastq"
simulated_sam_name = "_art.sam"

simulated_fq_1_name = "_art1.fastq"
simulated_fq_2_name = "_art2.fastq"

# K65R, K103N
pcr_error_positions = [650, 774]

package_dir = seq2simulate.__path__[0]
contamination_dir = os.path.join(package_dir, 'contamination')
human_file = os.path.join(contamination_dir, 'human_chr21')
env_file = os.path.join(contamination_dir, 'hxb2_env')

def merge_files_single(filenames, out_filename):
    """
    Merge the contents of a list of files into a single file

    Args:
        filenames: The filenames to merge.
        out_filename: Where to put the merged file.

    """

    with open(out_filename, 'w') as out_handle:
        for filename in filenames:
            with open(filename, 'r') as in_handle:
                for line in in_handle:
                    out_handle.write(line)

def merge_files(filenames, out_filename, paired_end):
    """
    Merge simulated files into one

    Args:
        filenames: The filenames to merge.
        out_filename: Where to put the merged file.
        paired_end: do we have pairs of files?

    """
    if paired_end:
        merge_files_single([f[0] for f in filenames], out_filename[0])
        merge_files_single([f[1] for f in filenames], out_filename[1])
    else:
        merge_files_single(filenames, out_filename)


def merge_sam_files(filenames, out_filename):
    """
    Merge the contents of a list of sam files into a single file, using the
    header from the first file.

    Args:
        filenames: The filenames to merge.
        out_filename: Where to put the merged file.

    """

    sq = []
    for filename in filenames:
        with open(filename, 'r') as handle:
            for line in handle:
                if line.startswith('@SQ'):
                    sq.append(line)
                if not line.startswith('@'):
                    break

    with open(out_filename, 'w') as out_handle:

        with open(filenames[0], 'r') as in_handle:
                out_handle.write(in_handle.readline())

        for header in sq:
            out_handle.write(header)

        for filename in filenames:
            with open(filename, 'r') as in_handle:
                for line in in_handle:
                    # ignore header lines from the other files
                    if line[0] != '@':
                        out_handle.write(line)

def split_into_two_sams(sam_filename):
    sam_1_filename = sam_filename + ".1.sam"
    sam_2_filename = sam_filename + ".2.sam"
    with open(sam_filename, 'r') as sam_original, \
        open(sam_1_filename, 'w') as sam_1, \
        open(sam_2_filename, 'w') as sam_2:
        while True:
            line1 = sam_original.readline()
            line2 = sam_original.readline()
            if not line2:
                break

            sam_1.write(line1)
            sam_2.write(line2)
    return sam_1_filename, sam_2_filename

def copy_sam_header(original_filename, new_handle):
    """
    Copy a SAM header from one file to another.

    Args:
        original_filename: The name of the original SAM file.
        new_handle: The handle to write to.

    """
    original_handle = open(
        original_filename, 'r')

    header = []
    next_header = original_handle.next()
    while next_header.startswith('@'):
        new_handle.write(next_header)
        next_header = original_handle.next()
    original_handle.close()


def unaligned_sam_file(original_sam_filename, working_dir):
    """
    Remove the alignment from a SAM format file.

    Args:
        original_sam_file: The aligned SAM file.
        working_dir: Put temp files here.

    """

    corrected_sam_filename = os.path.join(working_dir, str(uuid.uuid4()))
    corrected_sam_file = open(corrected_sam_filename, 'w')

    copy_sam_header(original_sam_filename, corrected_sam_file)

    incorrect_sam_file = custom_generators.parse_sam(
        open(original_sam_filename, 'r')
    )

    for read in incorrect_sam_file:
        
        read[custom_generators.CIGAR] = ''
        # 4 is the SAM flag signal for "unmapped"
        read[custom_generators.FLAG] = '4'
        corrected_sam_file.write('\t'.join(read))

    incorrect_sam_file.close()
    corrected_sam_file.close()

    shutil.move(corrected_sam_filename, original_sam_filename)


def correct_sam_file(original_sam_file, working_dir):
    """
    Replace the '=' and 'X' SAM arguments with 'M', because pysam doesn't
    understand '=' and 'X'.

    Args:
        original_sam_file: The file to change.
        working_dir: Put temp files here.
    """

    corrected_sam_filename = os.path.join(working_dir, str(uuid.uuid4()))
    corrected_sam_file = open(corrected_sam_filename, 'w')

    copy_sam_header(original_sam_file, corrected_sam_file)

    incorrect_sam_file = custom_generators.parse_sam(
        open(original_sam_file, 'r')
    )

    for read in incorrect_sam_file:
        if read[custom_generators.CIGAR] is not None \
            and read[custom_generators.CIGAR] != "":

            read[custom_generators.CIGAR] = \
                read[custom_generators.CIGAR].replace("=", "M")
            read[custom_generators.CIGAR] = \
                read[custom_generators.CIGAR].replace("X", "M")

        corrected_sam_file.write('\t'.join(read))

    incorrect_sam_file.close()
    corrected_sam_file.close()

    shutil.move(corrected_sam_filename, original_sam_file)


def run_art(
    sequence, platform, coverage, paired_end, working_dir, 
    forward_only=False, reverse_only=False
):
    """
    A wrapper for running the ART simulator from a SeqRecord object.

    Args:
        sequence_: The SeqRecord to simulate.
        platform: A Platform object.
        coverage: The required coverage.
        working_dir: Where to place temporary files.

    Returns:
        The filename of the simulated data
        
    """

    temp_file = os.path.join(working_dir, str(uuid.uuid4()))
    with open(temp_file, 'w') as out_handle:
        Bio.SeqIO.write([sequence], out_handle, 'fasta')

    result_file, sam_file = art.simulate(
                                temp_file, platform, coverage,
                                paired_end, working_dir
                            )
    print '.',
    os.unlink(temp_file)
    
    
    return result_file, sam_file


def run_art_with_pcr_error(
    sequence, platform, coverage, paired_end, working_dir
):
    """
    A wrapper for running the ART simulator from a SeqRecord object, first
    adding a PCR error.

    Args:
        sequence: The SeqRecord to simulate.
        platform: A Platform object.
        coverage: The required coverage.
        working_dir: Where to place temporary files.

    Returns:
        The filename of the simulated data
        
    """

    forward_clean_fq, forward_clean_sam = run_art(
        sequence, platform, coverage/4, paired_end,
        working_dir, forward_only=True
    )

    reverse_clean_fq, reverse_clean_sam = run_art(
        sequence, platform, coverage/2, paired_end,
        working_dir, reverse_only=True
    )

    pcr_error_sequence = copy.deepcopy(sequence)

    for error_pos in pcr_error_positions:
        error_to_add = random.choice(['A','C','G','T'])
        seq_list = list(str(pcr_error_sequence.seq))
        seq_list.insert(error_pos, error_to_add)
        pcr_error_sequence.seq = Bio.Seq.Seq(
            ''.join(seq_list), 
            Bio.Alphabet.Alphabet()
        )

    forward_dirty_fq, forward_dirty_sam = run_art(
        pcr_error_sequence, platform, coverage/4, 
        paired_end, working_dir, forward_only=True
    )

    if paired_end:
        out_file = [
            os.path.join(working_dir, str(uuid.uuid4())),
            os.path.join(working_dir, str(uuid.uuid4()))
        ]
        
    else:
        out_file = os.path.join(working_dir, str(uuid.uuid4()))
        

    merge_files(
        [forward_clean_fq, reverse_clean_fq, forward_dirty_fq], out_file,
        paired_end
    )

    out_sam = os.path.join(working_dir, str(uuid.uuid4()))
    merge_sam_files(
        [forward_clean_sam, reverse_clean_sam, forward_dirty_sam], out_sam
    )

    return out_file, out_sam


def simulate(sequence, sequence_file, platform, paired_end, working_dir):
    """
    Simulate reads from a sequencing platform from a multi-fasta file 
    containing the original sequence population to simulate.

    Args:
        sequence_file: The file from which to simulate.
        platform: A Platform object.
        working_dir: Where to place temporary files.

    Returns:
        The filename of the final simulated error file.
        
    """

    fq_filenames = []
    sam_filenames = []

    # generate twice as much coverage as required
    with open(sequence_file, 'r') as handle:
        # this is a tiny file, so it's ok to put it in memory
        records = list(Bio.SeqIO.parse(handle, 'fasta'))
        coverage = platform.coverage/len(records)

        coverage_multiplier = 4
        if platform == plat.roche:
            coverage_multiplier = 32

        print "Building", coverage_multiplier, "times coverage."

        for i, record in enumerate(records):
            if sequence.pcr_error:
                fq_file, sam_file = run_art_with_pcr_error(
                    record, platform, coverage * coverage_multiplier, 
                    paired_end, working_dir)
            else:
                fq_file, sam_file = run_art(
                    record, platform, coverage * coverage_multiplier, 
                    paired_end, working_dir)

            fq_filenames.append(fq_file)
            sam_filenames.append(sam_file)

    errors = [
        (sequence.human_error, human_file),
        (sequence.env_error, env_file)
    ]
    
    for error, error_file in errors:
        # simulate some reads from a human chromosome and add them
        if error:
            contaminate = random.uniform(0.1, 0.5)
            simulated_fq, simulated_sam = art.simulate(
                error_file, platform, int(platform.coverage*contaminate), 
                paired_end, working_dir
            )
            unaligned_sam_file(simulated_sam, working_dir)
            fq_filenames.append(simulated_fq)
            sam_filenames.append(simulated_sam)

    if paired_end:
        out_fq = (sequence_file + simulated_fq_1_name, 
            sequence_file + simulated_fq_2_name)
    else:
        out_fq = sequence_file + simulated_fq_name

    out_sam = sequence_file + simulated_sam_name
    merge_files(fq_filenames, out_fq, paired_end)
    merge_sam_files(sam_filenames, out_sam)
    correct_sam_file(out_sam, working_dir)

    if paired_end:
        out_sam = split_into_two_sams(out_sam)

    for filename in sam_filenames:
        os.unlink(filename)

    if(paired_end):
        for filenames in fq_filenames:
            os.unlink(filenames[0])
            os.unlink(filenames[1])
    else:
        for filename in fq_filenames:
            os.unlink(filename)

    print
    return out_fq, out_sam


def shuffle_and_assign_error_types(sequences, pcr_error=True):
    """
    Shuffle sequences, copy them, and assign certain common sequencing errors 
    to the shuffled sequences 0-6.  If less than 3 sequences are passed, no 
    errors will be assigned.  If less than 6 are passed, only one sequence 
    will be assigned each error.  Otherwise two sequences will be assigned 
    each of the error types.

    Args:
        sequences: the list of sequences to shuffle.    
    """
    
    random.shuffle(sequences)

    if len(sequences) >= 3 and len(sequences) < 6:
        if pcr_error:
            sequences[0].pcr_error = True
        sequences[1].env_error = True
        sequences[2].human_error = True
    elif len(sequences) >= 6:
        if pcr_error:
            sequences[3].pcr_error = True
        sequences[4].env_error = True
        sequences[5].human_error = True
