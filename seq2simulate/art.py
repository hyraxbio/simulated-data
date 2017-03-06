import os
import shutil
import subprocess
import uuid
import re

import Bio

import platform as plat
from platform import Platform
from maf2sam import maf2sam
import seq2simulate

def get_art_version(art_executable):
    c = subprocess.Popen(art_executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    art_text, error = c.communicate()
    return [int(i) for i in re.findall(r'Version (\d+\.\d+\.\d+)', art_text)[0].split('.')]

def sequence_length(sequence_file):
    """ Open a file and return the length of the first sequence.

        Args:
            sequence_file: Filename of the sequence.
    """
    with open(sequence_file, 'r') as handle:
        return len(list(Bio.SeqIO.parse(handle, "fasta"))[0])

def simulate(sequence_file, platform, coverage, paired_end, working_dir, extra_art_args=None):

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

    extra_args = {}
    if extra_art_args is not None:
        extra_args_separated = [str(i) for i in extra_art_args.split(',') ]
        pars = extra_args_separated[::2]
        vals = extra_args_separated[1::2]
        for par, val in zip(pars, vals):
            extra_args[par] = val 
 
    if platform == plat.illumina:
        # catch for version-specific params which changed in GreatSmokeyMountains
        art_version = get_art_version('art_illumina')
        ms = 'MS'
        if art_version[1] > 3:
            ms = 'MSv3'
        args = [
            'art_illumina', 
            '-sam',
            '-na', 
            '-i', sequence_file,
            '-l', str(platform.mean_read_length),
            # we simulate illumina miseq reads 
            '-ss', ms, 
            '-f', str(coverage),
            '-o', out_file,
        ]

        if paired_end:
            args.extend([
                '-p',
                '-m', '400',
                '-s', '50'
            ])

    elif platform == plat.pacbio_ccs:
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
    elif platform == plat.pacbio_clr:
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

    for par, val in extra_args.iteritems():
        args.append(par)
        args.append(val)
            
    devnull = open(os.devnull, 'w')
    subprocess.check_call(args, stdout=devnull, stderr=devnull)
    devnull.close()

    if platform == plat.pacbio_ccs or platform == plat.pacbio_clr:
        shutil.move(out_file + '_0001.fastq', out_file + '.fq')
        shutil.move(out_file + '_0001.ref', out_file + '.ref')
        shutil.move(out_file + '_0001.maf', out_file + '.maf')
        maf2sam(out_file + '.maf', out_file + '.ref')
        os.unlink(out_file + '.maf')
        os.unlink(out_file + '.ref')

    if platform == plat.roche or platform == plat.ion:
        os.unlink(out_file + '.stat')

    if paired_end:
        return (out_file + '1.fq', out_file + '2.fq'), out_file + '.sam'
    else:
        return out_file + '.fq', out_file + '.sam'
