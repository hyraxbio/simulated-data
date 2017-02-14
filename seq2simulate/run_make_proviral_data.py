from Bio import SeqIO, Seq, Alphabet
import os

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

def run_proviral(sequences, working_dir, out_dir, platform, paired_end, proviral_fraction):
    bio_sequences = [s for s in SeqIO.parse(sequences, 'fasta')]
    hypermutated_sequence_file = hypermutate_sequences(bio_sequences, working_dir)

    bio_sequences = [s for s in SeqIO.parse(sequences, 'fasta')]
    hyper_sequences = [s for s in SeqIO.parse(hypermutated_sequence_file, 'fasta')]
    # inject fraction of proviral sequences
    reduced_proviral_fraction = int(10*proviral_fraction) // 1

    bio_sequences *= 10-reduced_proviral_fraction
    hyper_sequences *= reduced_proviral_fraction
    total_sequences = bio_sequences + hyper_sequences

    full_filename = os.path.join(
        working_dir, 
        "mixed_hyperdata.fasta",
    )
    SeqIO.write(total_sequences, full_filename, 'fasta')

    print(len(total_sequences))
    platf = getattr(plat, platform)
    fastq_file, sam_file = art.simulate(full_filename, platf, platf.coverage, paired_end, out_dir)
    print 'Output saved in:', out_dir
