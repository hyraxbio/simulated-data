from Bio import SeqIO, Seq, Alphabet
import os
import random

from hypermutation import hypermutate
import platform as plat
import art

# proviral hypermutations per hundred bp
hypermutation_rate = 3

def hypermutate_sequences(sequences, working_dir):
    print('\n-------------------------------------------------------')
    print('Hypermutating evolved sequences (rate = {} per 100 bp).'.format(hypermutation_rate))
    print('-------------------------------------------------------\n')
    hypermutation_rates = [int(hypermutation_rate * len(str(s.seq)) // 100) for s in sequences]
    hyper_evolved_sequences = hypermutate.mutate_sequences([str(s.seq) for s in sequences], hypermutation_rates)
    for i, hseq in enumerate(hyper_evolved_sequences):
        sequences[i].seq = Seq.Seq(hseq, alphabet=Alphabet.SingleLetterAlphabet())
    full_filename = os.path.join(
        working_dir, 
        "hyperdata.fasta",
    )
    SeqIO.write(sequences, full_filename, 'fasta')
    return full_filename

def parse_fastq(filename):
    try:
        with open(filename, 'r') as f:
            fq = f.read()
        fq = fq.split('\n')
        fq = ['\n'.join(fq[i*4:i*4+4]) for i in range(len(fq)/4)]
        return fq
    except IOError:
        print('File not found.')


def run_proviral(sequences, working_dir, out_dir, platform, paired_end, proviral_fraction):
    bio_sequences = [s for s in SeqIO.parse(sequences, 'fasta')]
    hypermutated_sequence_file = hypermutate_sequences(bio_sequences, working_dir)

    platf = getattr(plat, platform)

    fastq_file0, sam_file0 = art.simulate(sequences, platf, platf.coverage, paired_end, out_dir)
    fastq_file1, sam_file1 = art.simulate(hypermutated_sequence_file, platf, platf.coverage, paired_end, out_dir)

    fq0 = parse_fastq(fastq_file0) 
    fq1 = parse_fastq(fastq_file1) 

    n_reads = min(len(fq0), len(fq1))
    n_proviral_reads = int(proviral_fraction * n_reads)

    print('Sampling {0} proviral reads out of {1} total reads (r={2:.3f}).'.format(n_proviral_reads, n_reads, n_proviral_reads/float(n_reads)))

    mixed_fastq = []
    while len(mixed_fastq) < n_proviral_reads:
        i = random.randint(0, len(fq1) - 1)
        mixed_fastq.append(fq1.pop(i))
    while len(mixed_fastq) < n_reads:
        i = random.randint(0, len(fq0) - 1)
        mixed_fastq.append(fq0.pop(i))
   
    mixed_fastq = '\n'.join(mixed_fastq)

    full_filename = os.path.join(
        out_dir, 
        "mixed_hyperdata.fastq",
    )
    with open(full_filename, 'w') as f:
        f.write(mixed_fastq) 

    print 'Output saved in:', out_dir
