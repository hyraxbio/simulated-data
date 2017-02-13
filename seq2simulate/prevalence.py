import itertools
import os
import random
import shutil
import uuid

import Bio
import pysam
import re

import art
import custom_generators
import drm
import sequencing_error

# regular expression that splits a cigar string into tuples
cigar_re = re.compile("([0-9]+)([A-Z])")

QUALITY_CUTOFF = 19
FASTQ_OFFSET = 33

ABOVE = 0
BELOW = 1

# start at RT 89
rt_no_coverage_start = ((89 + 155) * 3 - 3)
# end at RT 349
rt_no_coverage_end = ((349 + 155) * 3)

def range_resistant(platform, required_prevalence):
    """
    Considering the platform error rate, should we try to test this as 
    resistant?
    Args:
        platform: The sequencing Platform (class).
        prevalence: The DRM prevalence
    Returns:
        True if prevalence is well above error range, False otherwise
    """
    if required_prevalence >= platform.prevalence_error * 2:
        return True
    else:
        return False

def range_susceptible(platform, required_prevalence):
    """
    Considering the platform error rate, should we try to test this as 
    susceptible?
    Args:
        platform: The sequencing Platform (class).
        prevalence: The DRM prevalence.
    Returns:
        True if prevalence is well below error range, False otherwise.
    """

    if required_prevalence <= platform.prevalence_error / 2:
        return True
    else:
        return False

def range_unusable(platform, required_prevalence):
    """
    Is this prevalence in the no-man's land where no prevalence cutoff is
    guaranteed to produce either a susceptible or a resistant result?
    Args:
        platform: The sequencing Platform (class).
        prevalence: The DRM prevalence.
    Returns:
        True if prevalence is very close to error rate, meaning we cannot
        confidently say we'll call at resistant or susceptible.
    """
    return not range_resistant(platform, required_prevalence) \
           and not range_susceptible(platform, required_prevalence)

def sam_interval(sam):
    """
    What interval does a cigar string truly cover?

    Args:
        sam: The SAM entry.
    Returns:
        A tuple of the true start and end positions of the read.

    """

    if sam[custom_generators.FLAG] == 4:
        return False

    pos = int(sam[custom_generators.POS])
    length = len(sam[custom_generators.SEQ])

    cigar_tuples = cigar_re.findall(sam[custom_generators.CIGAR])
    ref_pos = pos
    query_pos = 0
    final_query_pos = -1

    for num_str, kind in cigar_tuples:

        num = int(num_str)
        if kind == 'I':
            length -= num
            query_pos += num
        elif kind == 'D':
            length += num
            ref_pos += num
        else:
            query_pos += num
            ref_pos += num

    return (pos, pos + length)

def drm_in_interval(drm, interval):
    """
    Is the drm in this interval?

    Args:
        drm: the drm in question
        interval [start, end] (inclusive)
    Returns:
        True or False

    """
    start, end = interval
    return drm.nucleotide_pos >= start and drm.nucleotide_pos <= end

def in_removed_interval(read):
    """
    Does the sam file stick out into the interval we are removing?
    """

    ret_val = False
    start, end = sam_interval(read[1])
    if end >= rt_no_coverage_start and start <= rt_no_coverage_end:
        ret_val = True

    if len(read) == 4:
        start, end = sam_interval(read[3])
        if end >= rt_no_coverage_start and start <= rt_no_coverage_end:
            ret_val &= True

    return ret_val

def can_add(sam, drm_dict, reads_per_drm, remove_rt):
    """
    A helper function for within_required_prevalence.  Are we allowed to
    add this DRM?

    Args:
        sam: A sam interval
        drm_dict: the drm dictionary
        remove_rt: are we removing a part of the reverse transcriptase gene?

    Returns:
        True/False
    """

    interval = sam_interval(sam)

    if remove_rt and in_removed_interval([0, sam]):
        return False

    # see if we can add this read without going over the number needed
    # for prevalence of any DRM.
    for drm, count in drm_dict.iteritems():
        if drm_in_interval(drm, interval) and count >= reads_per_drm:
                return False
    return True

def within_required_prevalence(read, drm_dict, reads_per_drm,
    remove_rt = False):
    """
    Checks whether a read should be added based on the values in a dictionary
    of DRMs.  Numbers of reads containing DRMs are added precisely to ensure
    counts are correct.  Other reads are added on a sliding probability scale
    to make up approximately correct coverage.

    Args:
        read: A read, sam pair
        drm_dict: A dictionary of DRMs of interest and the number of reads
                  added so far that contain each of those DRMs.
                  NB NB: drm_dict IS ALSO UPDATED BY THIS FUNCTION.
        reads_per_drm: The exact number of reads required per DRM.

    Returns:
        True if the read should be added, False otherwise.
    """

    # don't add any more reads if we've reached coverage for all drms.
    if not any([val < reads_per_drm for val in drm_dict.values()]):
        return False

    add_read = can_add(read[1], drm_dict, reads_per_drm, remove_rt)
    # indicates we've zipped together two files for paired end
    if len(read) == 4:
        add_read &= can_add(read[3], drm_dict, reads_per_drm, remove_rt)
    
    if not add_read:
        return False

    interval = sam_interval(read[1])
    # if we don't violate counts, go round again and increment counts
    contains_drm = False
    for drm, count in drm_dict.iteritems():
        if drm_in_interval(drm, interval):
                drm_dict[drm] += 1
                contains_drm = True

    if len(read) == 4:
        interval2 = sam_interval(read[3])
        for drm, count in drm_dict.iteritems():
            if (
                drm_in_interval(drm, interval2) and 
                    not drm_in_interval(drm, interval)
            ):
                drm_dict[drm] += 1
                contains_drm = True

    if contains_drm == True:
        return True
        
    else:

        return not any([val == reads_per_drm for val in drm_dict.values()])


def drms_in_read(drm_list, read):
    """
    Return which of a list of DRMs a read contains.

    Args:
        drm_list: a list of DRM objects.
        sam: a sam object

    Returns:
        A list of DRM objects contained by the read.
    """

    interval = sam_interval(read["sam"])
    result = [d for d in drm_list if drm_in_interval(d, interval)]

    # indicates a paired-end read: two fastqs and two sams
    if len(read) == 4:
        interval2 = sam_interval(read["sam2"])
        result2 = [d for d in drm_list if drm_in_interval(d, interval2)]
        result = list(set(result + result2))

    return result

class Blocker:
    """
    A little class containing information about reads that must be deleted
    in order to fulfill prevalence promises.
    """

    def __init__(self, drm):
        self.drm = drm
        self.blocking_drms = []
        self.num_to_delete = 0
        self.num_deleted = 0
        self.num_added = 0

def single_generator(fq_handle, sam_handle, drm_list, 
    reads_per_drm, remove_rt):

    """
    Return a generator to an output file and an associated dictionary
    to be filled.

    Args:
        fq_handle: an open file handle to a fastq file.
        sam_handle: an open file handle to a sam file.
        drm_list: a list of DRMs that must have exact prevalences.
        reads_per_drm: the number of reads to add per DRM.
        remove_rt: are we removing (almost) all RT DRMs?
    Returns:
        a tuple of:
            * an unpopulated dictionary of {DRM: coverage count}
            * the generator which, once exhausted, populates the dictionary
              and yields the reads to include.
    """

    drm_dict = { d : 0 for d in drm_list }

    in_generator = itertools.izip(
        custom_generators.parse_fastq(fq_handle), 
        custom_generators.parse_sam(sam_handle)
    )
    out_generator = ({
            "read": read[0], 
            "sam": read[1]
        } for read in in_generator \
            if within_required_prevalence(
                read, drm_dict, reads_per_drm, remove_rt
            ))
    return drm_dict, out_generator

def paired_generator(fq_handle1, fq_handle2, sam_handle1, sam_handle2, 
    drm_list, reads_per_drm, remove_rt):

    """
    Return a generator to an output file and an associated dictionary
    to be filled.

    Args:
        fq_handle1: first open file handle to fastq file.
        fq_handle2: second open file handle to fastq file.
        sam_handle1: first open file handle to sam file.
        sam_handle2: second open file handle to sam file.
        drm_list: a list of DRMs that must have exact prevalences.
        reads_per_drm: the number of reads to add per DRM.
        remove_rt: are we removing (almost) all RT DRMs?
    Returns:
        a tuple of:
            * an unpopulated dictionary of {DRM: coverage count}
            * the generator which, once exhausted, populates the dictionary
              and yields the reads to include.
    """

    drm_dict = { d : 0 for d in drm_list }

    in_generator = itertools.izip(
        custom_generators.parse_fastq(fq_handle1),
        custom_generators.parse_sam(sam_handle1),
        custom_generators.parse_fastq(fq_handle2),
        custom_generators.parse_sam(sam_handle2)
    )
    out_generator = ({
            "read": read[0],
            "sam": read[1],
            "read2": read[2],
            "sam2": read[3]
        } for read in in_generator \
            if within_required_prevalence(
                read, drm_dict, reads_per_drm, remove_rt
            ))
    return drm_dict, out_generator

def generator(fq_filename, sam_filename, drm_list, 
    reads_per_drm, remove_rt, paired_end):
    """
    Return a generator to an output file and an associated dictionary
    to be filled, as well as a set of open file handles

    Args:
        fq_handle: an open file handle to a fastq file.
        sam_handle: an open file handle to a sam file.
        drm_list: a list of DRMs that must have exact prevalences.
        reads_per_drm: the number of reads to add per DRM.
        remove_rt: are we removing (almost) all RT DRMs?
        paired_end: is this paired end data?
    Returns:
        a tuple of:
            * file handles to close
            * an unpopulated dictionary of {DRM: coverage count}
            * the generator which, once exhausted, populates the dictionary
              and yields the reads to include.
    """
    handles_to_close = []

    if paired_end:
        fq_handle1 = open(fq_filename[0], 'r')
        fq_handle2 = open(fq_filename[1], 'r')
        sam_handle1 = open(sam_filename[0], 'r')
        sam_handle2 = open(sam_filename[1], 'r')
        drm_dict, out_generator = paired_generator(
            fq_handle1, fq_handle2, sam_handle1, sam_handle2, drm_list, 
            reads_per_drm, remove_rt
        )
        handles_to_close.append(fq_handle1)
        handles_to_close.append(fq_handle2)
        handles_to_close.append(sam_handle1)
        handles_to_close.append(sam_handle2)
    else:
        fq_handle = open(fq_filename, 'r')
        sam_handle = open(sam_filename, 'r')
        drm_dict, out_generator = single_generator(
            fq_handle, sam_handle, drm_list, reads_per_drm, remove_rt
        )

        handles_to_close.append(fq_handle)
        handles_to_close.append(sam_handle)
    return handles_to_close, drm_dict, out_generator

def delete_generator(old_generator, drm_list, all_blockers):

    """
    Take in an unexhausted generator and return a new generator which
    removes reads which block fulfillment of DRM prevalence requirements.

    Args:
        old_generator: the original generator.
        drm_list: a list of DRMs that must have exact prevalences.
        all_blockers: a list of Blocker objects, which record how many
                      reads must be deleted containing each blocking drm set.
    Returns:
        yields reads, no return type.
    """

    for read in old_generator:
        do_yield = True
        for b in [b for b in all_blockers
            if b.num_to_delete > b.num_deleted]:
                
            if (set(drms_in_read(drm_list, read))
                == set(b.blocking_drms)):
                do_yield = False
                b.num_deleted += 1
                
        if do_yield:
            yield read

def add_generator_internal(in_generator, drm_list, blocker, remove_rt):

    """
    Take in a generator that zips a fastq and sam file and return a set
    of reads to add that replace reads previously deleted because they block
    fulfillment.
    Used by add_generator exclusively.

    Args:
        in_generator: a clean generator of fastq, sam pairs.
        drm_list: a list of DRMs that must have exact prevalences.
        blocker: the current Blocker, with populated num_deleted (otherwise
                 this generator does nothing)
        remove_rt: are we removing (almost) all RT DRMs?
    Returns:
        yields reads, no return type.
    """

    b_drms = blocker.blocking_drms
    b_drms.append(blocker.drm)
    for read in in_generator:
        if blocker.num_added == blocker.num_deleted:
            break

        if len(read) == 4:
            read_dict = {"read": read[0], "sam": read[1],
                "read2": read[2], "sam2": read[3]}
        else:
            read_dict = {"read": read[0], "sam": read[1]}
        if set(drms_in_read(drm_list, read_dict)) == \
            set(b_drms) and not \
            (remove_rt and in_removed_interval(read)):
                yield read_dict
                blocker.num_added += 1


def add_generator(fq_filename, sam_filename, drm_list, all_blockers, 
    remove_rt, paired_end):

    """
    Refulfill all prevalences by adding reads.  Designed to be chained to
    a delete_generator.

    Args:
        fq_filename: a fastq filename.
        sam_filename: a sam filename.
        drm_list: a list of DRMs that must have exact prevalences.
        all_blockers: blocking DRM sets, with populated num_deleted.
        remove_rt: are we removing (almost) all RT DRMs?
    Returns:
        A tuple containing:
            * A generator that adds reads to re-fulfill all prevalences.
            * A list of open file handles (NB because this generator is not
              exhausted by this function, file handles must be closed
              only when generator is exhausted).
    """

    final_generator = iter([])
    handles_to_close = []

    for blocker in all_blockers:
        if paired_end:
            fq_handle1 = open(fq_filename[0], 'r')
            fq_handle2 = open(fq_filename[1], 'r')
            sam_handle1 = open(sam_filename[0], 'r')
            sam_handle2 = open(sam_filename[1], 'r')
            in_generator = itertools.izip(
                custom_generators.parse_fastq(fq_handle1),
                custom_generators.parse_sam(sam_handle1),
                custom_generators.parse_fastq(fq_handle2),
                custom_generators.parse_sam(sam_handle2)
            )
        else:
            fq_handle = open(fq_filename, 'r')
            sam_handle = open(sam_filename, 'r')
            handles_to_close.extend([fq_handle, sam_handle])
            in_generator = itertools.izip(
                custom_generators.parse_fastq(fq_handle), 
                custom_generators.parse_sam(sam_handle)
            )
        generator_frag = add_generator_internal(in_generator, drm_list,
            blocker, remove_rt)
        final_generator = itertools.chain(final_generator, generator_frag)

    return final_generator, handles_to_close

def next_blocking_drm(drm_dict, all_blockers, unfulfilled_drm):
    """
    Figure out which fulfilled DRM is currently blocking fulfillment of
    the unfulfilled DRM.  Sometimes more than one DRM is blocking fulfillment
    concurrently.

    Args:
        drm_dict: a fulfilled drm_dict that shows which DRMs don't have
                  complete prevalence yet.
        all_blockers: blocking DRMs already found.
        unfulfilled_drm: the DRM being blocked.
        
    Returns:
        The next DRM to add to the blocker list.
    """

    return min(
        [
            d for d, val in drm_dict.iteritems() if \
                d not in (
[drm for b in all_blockers for drm in (b.blocking_drms)] + [unfulfilled_drm]
        )],
        # min is defined as the next-closest distance to current drm 
        key=(
            lambda x : abs(
                x.nucleotide_pos - unfulfilled_drm.nucleotide_pos
            )
        )
    )

def blockers(fq_filename, sam_filename, drm_dict, 
            reads_per_drm, remove_rt, paired_end):

    """
    Refulfill all prevalences by adding reads.  Designed to be chained to
    a delete_generator.

    Args:
        fq_filename: a fastq filename.
        sam_filename: a sam filename.
        drm_dict: a previously-populated DRM dictionary.
        reads_per_drm: the number of reads to add per DRM.
        remove_rt: are we removing (almost) all RT DRMs?
        paired_end: paired end data?
    Returns:
        A set of DRM sets blocking fulfillment, and how many reads containing
        exactly that DRM set must be deleted to complete fulfillment.
    """

    all_blockers = []
    while any([val < reads_per_drm for val in drm_dict.values()]):

        # find the next unfulfilled drm prevalence
        unfulfilled_drm = [d for d, val in drm_dict.iteritems()
            if val < reads_per_drm][0]

        # how much prevalence is missing?
        prevalence_diff = reads_per_drm - drm_dict[unfulfilled_drm]
        
        current_blockers = {ABOVE: [], BELOW: []}

        while (sum([b.num_to_delete for b in (
            current_blockers[ABOVE] + current_blockers[BELOW])]) \
        ) < prevalence_diff:

            blocking_drm = next_blocking_drm(drm_dict, current_blockers[ABOVE]
                + current_blockers[BELOW], unfulfilled_drm)

            blocker = Blocker(unfulfilled_drm)
            position = ABOVE
            if blocking_drm.nucleotide_pos < unfulfilled_drm.nucleotide_pos:
                position = BELOW

            # we need all currently blocking drms to unblock
            if current_blockers[position]:
                blocker.blocking_drms = list(current_blockers[position][-1].blocking_drms)
            blocker.blocking_drms.append(blocking_drm)
            current_blockers[position].append(blocker)

            print "Fulfillment of", unfulfilled_drm, "is blocked by", \
                blocking_drm, \
                "- correcting a difference of", prevalence_diff
            print drm_dict

            handles_to_close, _, temp_generator = generator(fq_filename, 
                sam_filename, drm_dict.keys(), reads_per_drm, remove_rt,
                paired_end)

            for read in temp_generator:
                if sum([b.num_to_delete for b in (
                    current_blockers[ABOVE] + current_blockers[BELOW])]
                ) == prevalence_diff:
                    break

                for b in (current_blockers[ABOVE] + current_blockers[BELOW]):

                    if (set(drms_in_read(drm_dict.keys(), read)) 
                        == set(b.blocking_drms)):
                        b.num_to_delete += 1

            for handle in handles_to_close:
                handle.close()

        all_blockers.extend(current_blockers[ABOVE])
        all_blockers.extend(current_blockers[BELOW])

        for b in (current_blockers[ABOVE] + current_blockers[BELOW]):
            drm_dict[b.drm] += b.num_to_delete

    return all_blockers

def populate_dictionary(fq_filename, sam_filename, drm_list, 
        reads_per_drm, remove_rt, paired_end):

    """
    Make a first, pre-corrected pass at fulfilling all prevalence
    requirements.  Just returns the dictionary, so we can see how wrong
    we were.

    Args:
        fq_filename: fastq filename(s) - a tuple if we're in paired-end mode.
        sam_filename: a sam filename.
        drm_list: a list of DRMs that must have exact prevalences.
        reads_per_drm: the number of reads to add per DRM.
        remove_rt: are we removing (almost) all RT DRMs?
    Returns:
        The populated dictionary.
    """

    handles_to_close, drm_dict, out_generator = generator(
            fq_filename, sam_filename, drm_list, reads_per_drm, 
            remove_rt, paired_end
        )

    # populate drm_dict by running the generator until the dictionary
    # is full or there are no more reads
    for _ in out_generator:
        if not any([val < reads_per_drm for val in drm_dict.values()]):
            break

    for h in handles_to_close:
        h.close()

    return drm_dict

def select_sequences(
    out_file, fq_filename, sam_filename,
    drm_list, reads_per_drm, num_other,
    remove_rt, paired_end
):
    """
    Output reads unless they cause an incorrect prevalence.

    Args:
        out_file: The filename to output the sequences to.
        in_file: The raw file to read from.
        drm_list: A list of DRM objects of interest.
        reads_per_drm: Prevalence at drm positions.
        num_other: Total number of other nucleotides to add.
    """

    # Don't require fulfillment of those DRMs at positions
    # where we are removing coverage
    if remove_rt:
        drm_list = [d for d in drm_list if \
            not drm_in_interval(d, (
                rt_no_coverage_start, rt_no_coverage_end))
        ]

    # make a first pass at fulfilling all prevalences by running through
    # a generator which selects reads, until selecting more reads would
    # produce too much coverage at a position.
    drm_dict = populate_dictionary(fq_filename, sam_filename, drm_list, 
        reads_per_drm, remove_rt, paired_end)
    
    # the same generator as above, except this time we don't actually
    # run the generator: we're going to wrap it in two more generators.
    handles_to_close, _, out_generator = generator(fq_filename, sam_filename,
        drm_list, reads_per_drm, remove_rt, paired_end)


    # if the first pass at running through the generator doesn't produce
    # exactly the correct coverage at all DRM positions:
    if any([val < reads_per_drm for val in drm_dict.values()]):
        # figure out which DRMs are "blocking" the fulfillment of other
        # DRMs (usually the adjacent ones)
        all_blockers = blockers(fq_filename, sam_filename, drm_dict, 
            reads_per_drm, remove_rt, paired_end)
        # remove reads that contain *only* the blockers until the blocking DRM
        # and the unfulfilled DRM have the same coverage
        deleted_generator = delete_generator(out_generator, 
            drm_list, all_blockers)
        # add reads that contain *both* the blocker and the unfulfilled DRM
        # until they both have the same coverage.
        reads_to_add, add_handles = add_generator(fq_filename, sam_filename, 
            drm_list, all_blockers, remove_rt, paired_end)
        handles_to_close.extend(add_handles)

        # chain the adding generator to the deletion generator
        out_generator = itertools.chain(deleted_generator, reads_to_add)

    # run the final generator, which may or may not be a combination
    # of an initial generator, a generator that doesn't yield certain
    # reads, and a final generator that adds additional reads.
    if paired_end:
        with open(out_file[0], 'w') as out_handle1, \
             open(out_file[1], 'w') as out_handle2:
            for r in out_generator:
                out_handle1.write("%s" % ''.join(r["read"]))
                out_handle2.write("%s" % ''.join(r["read2"]))
    else:
        with open(out_file, 'w') as out_handle:
            for r in out_generator:
                out_handle.write("%s" % ''.join(r["read"]))

    for h in handles_to_close:
        h.close()

def append_files(out_file, in_file):
    """
    Append one file to another.

    Args:
        out_file: the file to append to.
        in_file: the file to be appended.
    """

    with open(in_file, 'r') as handle, open(out_file, 'a') as out_handle:
        out_handle.write(handle.read())

def tag_names(filename, sequence_name, working_dir):
    """
    Change all sequence names to a particular name.

    Args:
        filename: The file to rename.
        sequence_name: The name to call the sequence
        working_dir: Put temp files here.
    """

    temp_file = os.path.join(working_dir, str(uuid.uuid4()))

    with open(filename, 'r') as handle, open(temp_file, 'w') as out_handle:
        for record in custom_generators.parse_fastq(handle):
            record = (
                ''.join(['@', sequence_name, '\n']), 
                record[1], 
                record[2], 
                record[3]
            )
            out_handle.write(''.join(record))

    shutil.move(temp_file, filename)


def randomize_names(filename, working_dir):
    """
    Write out sequences with randomized names.

    Args:
        filename: The file to randomize.
        working_dir: Put temp files here.
    """

    temp_file = os.path.join(working_dir, str(uuid.uuid4()))

    with open(filename, 'r') as handle, open(temp_file, 'w') as out_handle:
        for record in custom_generators.parse_fastq(handle):
            record = (
                ''.join(['@', str(uuid.uuid4()), '\n']), 
                record[1], 
                record[2], 
                record[3]
            )
            out_handle.write(''.join(record))

    shutil.move(temp_file, filename)

def randomize_names_paired(filename1, filename2, working_dir):
    """
    Write out paired end sequences with randomized names.

    Args:
        filename1: The first paired file to randomize.
        filename2: The second paired file to randomize.
        working_dir: Put temp files here.
    """

    temp_file1 = os.path.join(working_dir, str(uuid.uuid4()))
    temp_file2 = os.path.join(working_dir, str(uuid.uuid4()))

    with open(filename1, 'r') as handle1, \
            open(filename2, 'r') as handle2, \
            open(temp_file1, 'w') as out_handle1, \
            open(temp_file2, 'w') as out_handle2:
        final_generator = itertools.izip(
            custom_generators.parse_fastq(handle1),
            custom_generators.parse_fastq(handle2)
        )
        for record1, record2 in final_generator:
            name = str(uuid.uuid4())
            record1 = (
                ''.join(['@', name, "/1", '\n']), 
                record1[1], 
                record1[2], 
                record1[3]
            )
            record2 = (
                ''.join(['@', name, "/2", '\n']), 
                record2[1], 
                record2[2], 
                record2[3]
            )
            out_handle1.write(''.join(record1))
            out_handle2.write(''.join(record2))

    shutil.move(temp_file1, filename1)
    shutil.move(temp_file2, filename2)



def produce_prevalence(
                    prevalence, platform, paired_end, 
                    susceptible_file, susceptible_sam,
                    resistant_file, resistant_sam,
                    sequence, working_dir, out_dir
):
    """
    Takes in a file of susceptible fastqs and a file of resistant fastqs,
    and produces a file containing the correct prevalence of the
    resistant sequences at the drms in question.  Other sequences are added
    at random.

    Args:
        prevalence: The required prevalence at resistant positions. 
        platform: The sequencing platform from which reads were simulated.
        susceptible_file: The simulated susceptible read file.
        resistant_file: The simulated resistant read file.
        sequence: The Sequence object from which reads were simulated.
        out_dir: The final output folder.

    Returns:
        The name of the final file simulated.
    """

    print sequence.susceptible.id, '-', sequence.resistant.id, \
            'at', (prevalence*100), "% DRM prevalence"

    num_resistant = int(platform.coverage * prevalence)
    num_susceptible = int(platform.coverage * (1.0 - prevalence))

    num_total = (
        (len(sequence.susceptible) / platform.mean_read_length) \
        * platform.coverage
    )

    final_filename = os.path.join(out_dir, str(uuid.uuid4()))
    final_output_filename = final_filename
    if paired_end:
        final_filename = (
            final_filename + "_1.fastq", 
            final_filename + "_2.fastq"
        )
    else:
        final_filename += ".fastq"

    select_sequences(
        final_filename, susceptible_file, susceptible_sam, 
        sequence.drms, num_susceptible, 0, sequence.remove_rt, paired_end
    )

    final_resistant_filename = os.path.join(working_dir, str(uuid.uuid4()))
    if paired_end:
        final_resistant_filename = (
            final_resistant_filename + "_1", 
            final_resistant_filename + "_2"
        )

    select_sequences(
        final_resistant_filename, resistant_file, resistant_sam,
        sequence.drms, num_resistant, 0, sequence.remove_rt, paired_end
    )

    # add in some junk
    errors = [
        (sequence.human_error, sequencing_error.human_file),
        (sequence.env_error, sequencing_error.env_file)
    ]
    
    for error, error_file in errors:
        if error:
            print "Adding contamination."
            contaminate = random.uniform(0.1, 0.5)
            error_fq, _ = art.simulate(
                error_file, platform, int(float(platform.coverage)*contaminate), 
                paired_end, working_dir
            )
            if paired_end:
                append_files(final_filename[0], error_fq[0])
                append_files(final_filename[1], error_fq[1])
            else:
                append_files(final_filename, error_fq)

    if paired_end:
        append_files(final_filename[0], final_resistant_filename[0])
        os.unlink(final_resistant_filename[0])
        append_files(final_filename[1], final_resistant_filename[1])
        os.unlink(final_resistant_filename[1])
        randomize_names_paired(final_filename[0], final_filename[1],
            working_dir)
    else:
        # comment out sequence tags to "run blind":
        # tag_names(final_filename, "susceptible", working_dir)
        # tag_names(final_resistant_filename, "resistant", working_dir)
        append_files(final_filename, final_resistant_filename)
        os.unlink(final_resistant_filename)
        randomize_names(final_filename, working_dir)

    return final_output_filename



    
