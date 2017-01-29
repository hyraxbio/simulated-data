import copy
import itertools
import os
import random
import shutil
import subprocess
import sys
import uuid

import Bio
import pysam

import art

import custom_generators
import platform as plat
import drm as drm_mod
import seq2simulate

# values for ion torrent and miseq taken from
# http://aem.asm.org/content/early/2014/09/22/AEM.02206-14.full.pdf+html

simulated_fq_name = "_art.fastq"
simulated_sam_name = "_art.sam"

simulated_fq_1_name = "_art1.fastq"
simulated_fq_2_name = "_art2.fastq"

# approximate locations of K65R, K103N, taken from the end of integrase
# the gag-pol region is very variable, so better to count from end.
pcr_error_offsets = [2353, 2236]
# we'll only include pcr error when DRMs don't work.
pcr_error_inclusion_cutoff = 140

# flag numbers indicating forward and reverse
FORWARD_FLAG = 8
REVERSE_FLAG = 24

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
        handles = [open(f, 'r') for f in filenames]
        while True:
            all_eof = True
            for handle in handles:
                line = "".join([handle.readline() for _ in range(4)])
                if line is not None and line != "":
                    all_eof = False
                    out_handle.write(line)
            if all_eof:
                break
        for handle in handles:
            handle.close()

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

def select_only_flag(
    fastq_file, sam_file, working_dir, flag):
    """
    Select only fastq entries that have a particular SAM flag
    associated with them (usually forward or reverse).

    Args:
        fastq_file: the FASTQ file to take from
        sam_file: the SAM file to take from
        working_dir: the working directory
        flag: the flag to select

    Returns: (new fq file, new sam file)

    """

    out_fq_path = os.path.join(working_dir, str(uuid.uuid4()))
    out_sam_path = os.path.join(working_dir, str(uuid.uuid4()))

    with open(fastq_file, 'r') as fq_handle, \
         open(sam_file, 'r') as sam_handle:

        in_generator = itertools.izip(
            custom_generators.parse_fastq(fq_handle), 
            custom_generators.parse_sam(sam_handle)
        )

        in_generator, in_generator2 = itertools.tee(in_generator)

        with open(out_fq_path, 'w') as fq_out, \
             open(out_sam_path, 'w') as sam_out:

            for p in in_generator:
                if int(p[1][custom_generators.FLAG]) != flag:
                    continue
                fq_out.write("".join(p[0]))
                sam_out.write("\t".join(p[1]))

    return out_fq_path, out_sam_path

def run_art(
    sequence, platform, coverage, paired_end, working_dir, 
    forward_only=False, reverse_only=False
):
    """
    A wrapper for running the ART simulator from a SeqRecord object.

    Args:
        sequence: The SeqRecord to simulate.
        platform: A Platform object.
        coverage: The required coverage.
        working_dir: Where to place temporary files.
        forward_only: only simulate forward sequences.
        reverse_only: only simulate reverse sequences.

    Returns:
        The filename of the simulated data
        
    """

    assert not (forward_only and reverse_only)

    temp_file = os.path.join(working_dir, str(uuid.uuid4()))
    with open(temp_file, 'w') as out_handle:
        Bio.SeqIO.write([sequence], out_handle, 'fasta')

    result_file, sam_file = art.simulate(
                                temp_file, platform, coverage,
                                paired_end, working_dir
                            )
    os.unlink(temp_file)

    print '.',
    sys.stdout.flush()

    if forward_only:
        return select_only_flag(result_file, sam_file,
            working_dir, 0)
    elif reverse_only:
        return select_only_flag(result_file, sam_file,
            working_dir, 16)
    
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

    for error_pos in pcr_error_offsets:
        error_to_add = 'C'
        seq_list = list(str(pcr_error_sequence.seq))
        for i in range(-9, 9, 3):
            seq_list[len(seq_list) - error_pos + i] = error_to_add
        pcr_error_sequence = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(
            ''.join(seq_list), 
            Bio.Alphabet.Alphabet()
            ), id=sequence.id, name=sequence.name)

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
        [forward_clean_fq, forward_dirty_fq, reverse_clean_fq], out_file,
        paired_end
    )

    out_sam = os.path.join(working_dir, str(uuid.uuid4()))
    merge_sam_files(
        [forward_clean_sam, forward_dirty_sam, reverse_clean_sam], out_sam
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

    sequences[0].human_error = True
    sequences[1].env_error = True
    
    if pcr_error:
        pcr_added = False
        for sequence in sequences:
            can_add_pcr = True
            for drm in sequence.drms:
                if drm.locus == drm_mod.RT and not \
                ((drm.relative_pos > pcr_error_inclusion_cutoff) \
                or (drm.delete)):
                    can_add_pcr = False
                    break
            if can_add_pcr:
                sequence.pcr_error = True
                pcr_added = True
                break
        if not pcr_added:
            print "Failed to add PCR error."


