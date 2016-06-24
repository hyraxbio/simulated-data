import appdirs
import csv
import os
import json
import gzip
import shutil
import glob

from Bio import SeqIO, SeqRecord

def write_settings(folder, settings, samples):
    settings["sample_details"] = samples

    with open(os.path.join(folder, "settings.json"), 'w') as f:
        f.write(json.dump(settings, indent=2)) 



def compress(file):
    output = file + '.gz'
    with open(file, 'rb') as f_in, gzip.open(output, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    return output


def illumina_single(folder):
    """
    Package a folder of fastqs into a (roughly) Illumina Miseq compliant 
    folder (single-ended).

    Args:
        folder: The folder to package.

    Returns:
        True on completion.
    """
    if not os.path.isdir(folder):
        return

    settings = {"platform":1}
    samples = []
    for s in glob.glob(os.path.join(folder, "*.fastq")):
        file = os.basename(compress(s))
        sid = os.basename(s).replace('.fastq', '')
        samples.append({"sample_name":sid, "file_names":[file]})

    write_settings(folder, settings, samples)

    return True


def illumina_paired(folder):
    """
    Package a folder of fastqs into a (roughly) Illumina Miseq compliant 
    folder (paired-end).

    Args:
        folder: The folder to package.

    Returns:
        True on completion.
    """

    if not os.path.isdir(folder):
        return

    settings = {"platform":1, "paired":true}
    samples = []

    for s in glob.glob(os.path.join(folder, "*.fastq")):
        sid = os.basename(s).replace('.fastq', '')
        if sid.endswith('_1'):
            sid = sid.replace('_1', '')
            s2 = s.replace('_1.fastq', '_2.fastq')
            if os.path.isfile(s2):
                file1 = os.basename(compress(s))
                file2 = os.basename(compress(s2))
                samples.append({"sample_name":sid, "filenames":[file1, file2]})
            else:
                raise ValueError("Expecting pairs of files ending in _1 and _2")
    
    write_settings(folder, settings, samples)    
            
    return True


def ion(folder):
    """
    Package a folder of fastqs into an Ion Torrent compatible format

    Args:
        folder: The folder to package.

    Returns:
        True on completion.
    """

    if not os.path.isdir(folder):
        return

    settings = {"platform":0}
    samples = []
    for s in glob.glob(os.path.join(folder, "*.fastq")):
        file = os.basename(compress(s))
        sid = os.basename(s).replace('.fastq', '')
        samples.append({"sample_name":sid, "file_names":[file]})

    write_settings(folder, settings, samples)

    return True


roche_filename = "IULK3W101.fastq"
roche_details_preamble = "midname,patientname\n"

roche_prefix_characters = "tcag"

def read_generator(folder, file_list, mid_file):
    """
    Build combined reads with MIDs.

    Args:
        folder: The folder to package.
        file_list: A snapshot of the files in the folder at packaging start
        mids: A csv containing all allowed roche MIDs

    Returns:
        a tuple of sequence records and sample details
    """
    sequences = []
    samples = []
    
    with open(os.path.join(folder, roche_details_name), 'w') as details, \
        open(mid_file, 'r') as mid_handle:

        details.write(roche_details_preamble)
        mids = csv.reader(mid_handle)
        for seq_file, mid in zip(file_list, mids):
            in_seqs = SeqIO.parse(os.path.join(folder, seq_file), "fastq")
            prefix_string = roche_prefix_characters + mid[1]
            prefix_qual = [40] * 14
            sid = seq_file.replace('.fastq', '')
            samples.append({"sample_name":sid, "mid":mid[0]})
            for s in in_seqs:
                new_record = SeqRecord.SeqRecord(prefix_string + s.seq, 
                    id=s.id, 
                    name=s.name, 
                    description=s.description)
                new_record.letter_annotations["phred_quality"]=\
                    prefix_qual + s.letter_annotations["phred_quality"]
                sequences.append(new_record)
            os.unlink(os.path.join(folder, seq_file))

    return sequences, samples

def roche(folder, mid_file):
    """
    Package a folder of fastqs into a single Roche/454 file.

    Args:
        folder: The folder to package.
        mid_file: A csv containing all allowed roche MIDs

    Returns:
        True on completion.
    """

    settings = {"platform":2}
    
    if not os.path.isdir(folder):
        return

    file_list = os.listdir(folder)
    sequence, samples = read_generator(folder, file_list, mid_file)
    SeqIO.write(
        sequences, 
        os.path.join(folder, roche_filename), 
        "fastq"
    )
    write_settings(folder, settings, samples)

    return True


