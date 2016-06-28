import appdirs
import csv
import os
import gzip
import shutil
import glob

from Bio import SeqIO, SeqRecord

def compress(file):
    output = file + '.gz'
    with open(file, 'rb') as f_in, gzip.open(output, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(file)
    return output

illumina_details_name = "SampleSheet.csv"

illumina_details_preamble_single = """[Header]
IEMFileVersion,4
Investigator Name,JoeBloggs
Experiment Name,MY_HIV_EXPERIMENT
Date,2015/07/30
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,Nextera XT v2
Description,HIV amplicons 3.0KB
Chemistry,Amplicon

[Reads]
301

[Settings]
ReverseComplement,0
Adapter,CTGTCTCTTATACACATCT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
"""

illumina_details_preamble_paired = """[Header]
IEMFileVersion,4
Investigator Name,JoeBloggs
Experiment Name,MY_HIV_EXPERIMENT
Date,2015/07/30
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,Nextera XT v2
Description,HIV amplicons 3.0KB
Chemistry,Amplicon

[Reads]
301
301

[Settings]
ReverseComplement,0
Adapter,CTGTCTCTTATACACATCT

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
"""

illumina_details_postfix = "MY_HIV_270715,B09,N711-A,AAGAGGCA,S503-A,TATCCTCT,MY_HIVDR,"

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

    samples = []
    file_list = glob.glob(os.path.join(folder, "*.fastq"))
    with open(os.path.join(folder, illumina_details_name), 'w') as f:
        f.write(illumina_details_preamble_single)

        for s in file_list:
            pid = os.basename(s).replace('.fastq', '')
            f.write("%s,%s,%s\n" % (pid, pid, illumina_details_postfix))
            filename = os.basename(compress(path.join(folder, s)))
            samples.append({"sample_name": pid, "file_names": [filename]})

    return samples

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

    samples = []
    file_list = glob.glob(os.path.join(folder, "*.fastq"))
    with open(os.path.join(folder, illumina_details_name), 'w') as f:
        f.write(illumina_details_preamble_paired)
        for s in file_list:
            pid = os.basename(s).replace('.fastq', '')
            if pid.endswith('_1'):
                pid = pid.replace('_1', '')
                # only write out the first file of the pair
                f.write("%s,%s,%s\n" % (pid, pid, illumina_details_postfix))
                s2 = s.replace(s, "_1.fastq", "_2.fastq")
                if os.path.isfile(s2):
                    file1 = os.basename(compress(s))
                    file2 = os.basename(compress(s2))
                    samples.append({"sample_name":sid, "filenames":[file1, file2]})
                else:
                    raise ValueError("Expecting pairs of files ending in _1 and _2")
            
    return samples

ion_details_name = "PatientDetails.csv"
ion_details_preamble = "Sample ID,File Prefix\n"

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

    samples = []
    file_list = glob.glob(os.path.join(folder, "*.fastq"))
    with open(os.path.join(folder, ion_details_name), 'w') as f:
        f.write(ion_details_preamble)

        for s in file_list:
            pid = os.basename(s).replace('.fastq', '')
            f.write("%s,%s\n" % (pid, pid))
            filename = os.basename(compress(path.join(folder, s)))
            samples.append({"sample_name": pid, "file_names": [filename]})

    return samples


roche_filename = "IULK3W101.fastq"
roche_details_name = "MIDList.csv"
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
        True on completion.
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
            pid = seq_file.replace('.fastq', '')
            details.write("%s,%s\n" % (mid[0], pid))
            samples.append({"sample_name":pid, "mid":mid[0]})
            for s in in_seqs:
                new_record = SeqRecord.SeqRecord(prefix_string + s.seq, 
                    id=s.id, 
                    name=s.name, 
                    description=s.description)
                new_record.letter_annotations["phred_quality"]=\
                    prefix_qual + s.letter_annotations["phred_quality"]
                sequences.append(new_record)
            os.unlink(os.path.join(folder, seq_file))

    return sequences, sampes

def roche(folder, mid_file):
    """
    Package a folder of fastqs into a single Roche/454 file.

    Args:
        folder: The folder to package.
        mid_file: A csv containing all allowed roche MIDs

    Returns:
        True on completion.
    """
    
    if not os.path.isdir(folder):
        return

    file_list = os.listdir(folder)
    sequence, samples = read_generator(folder, file_list, mid_file)
    SeqIO.write(
        sequences, 
        os.path.join(folder, roche_filename), 
        "fastq"
    )
    return samples








