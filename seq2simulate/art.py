import os
import shutil
import subprocess
import uuid

import Bio

from platform import Platform
from maf2sam import maf2sam
import seq2simulate

package_dir = seq2simulate.__path__[0]
roche_profile_dir = os.path.join(package_dir, 'profiles/roche')
ion_profile_dir = os.path.join(package_dir, 'profiles/ion')
pacbio_ccs_profile = os.path.join(package_dir, 'profiles/model_qc_ccs')
pacbio_clr_profile = os.path.join(package_dir, 'profiles/model_qc_clr')

illumina = Platform(50000, 250)
roche = Platform(2000, 320, profile=roche_profile_dir)
ion = Platform(10000, 320, profile=ion_profile_dir)
pacbio_ccs = Platform(2000, 250, profile=pacbio_ccs_profile)
pacbio_clr = Platform(2000, 250, profile=pacbio_profile)

def sequence_length(sequence_file):
    """ Open a file and return the length of the first sequence.

        Args:
            sequence_file: Filename of the sequence.
    """
    with open(sequence_file, 'r') as handle:
        return len(list(Bio.SeqIO.parse(handle, "fasta"))[0])

def simulate(sequence_file, platform, coverage, paired_end, working_dir):

    """ Wrapper for the sequence simulators ART and pbsim.

        Args:
            sequence_file: Filename of the sequence from which to simulate
                        reads.
            platform: Platform object of the sequencing platform.
            working_dir: Where to put temporary files.

        Returns:
            (fastq_file, sam_file)
            The fastq_file and sam_file created.
    """
    out_file = os.path.join(working_dir, str(uuid.uuid4()))

    if platform == illumina:
        args = [
            'art_illumina', 
            '-sam',
            '-na', 
            '-i', sequence_file,
            '-l', str(platform.mean_read_length),
            # we simulate illumina miseq reads 
            '-ss', 'MS', 
            '-f', str(coverage), 
            '-o', out_file
        ]
        if paired_end:
            args.extend([
                '-p',
                '-m', '400',
                '-s', '50'
            ])

    elif platform == pacbio_ccs:
        args = [
            'pbsim',
            '--data-type',
            'CCS',
            '--model_qc',
            platform.profile,
            '--prefix',
            out_file, 
            '--depth',
            str(coverage),
            sequence_file,
        ]
    elif platform == pacbio_clr:
        args = [
            'pbsim',
            '--data-type',
            'CLR',
            '--model_qc',
            platform.profile,
            '--prefix',
            out_file, 
            '--depth',
            str(coverage),
            sequence_file,
        ]

    else:
        args = [
            'art_454',
            '-s',
            '-p', platform.profile,
            sequence_file,
            out_file,
            str(coverage)
        ]

    devnull = open(os.devnull, 'w')
    subprocess.check_call(args, stdout=devnull, stderr=devnull)
    devnull.close()

    if platform == pacbio:
        shutil.move(out_file + '_0001.fastq', out_file + '.fq')
        shutil.move(out_file + '_0001.ref', out_file + '.ref')
        shutil.move(out_file + '_0001.maf', out_file + '.maf')
        maf2sam(out_file + '.maf', out_file + '.ref')
        os.unlink(out_file + '.maf')
        os.unlink(out_file + '.ref')

    if platform == roche or platform == ion:
        os.unlink(out_file + '.stat')

    if paired_end:
        return (out_file + '1.fq', out_file + '2.fq'), out_file + '.sam'
    else:
        return out_file + '.fq', out_file + '.sam'
