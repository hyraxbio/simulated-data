import sys
from Bio import SeqIO

def split(sus, res):

    with open(sus, 'r') as sus_handle, open(res, 'r') as res_handle:

        for susceptible, resistant in zip(
            SeqIO.parse(sus_handle, "fasta"), 
            SeqIO.parse(res_handle, "fasta")
        ):
            with open(susceptible.id, 'w') as out_handle:
                SeqIO.write([susceptible], out_handle, "fasta")

            with open(resistant.id, 'w') as out_handle:
                SeqIO.write([resistant], out_handle, "fasta")


if __name__ == '__main__':

    split (sys.argv[1], sys.argv[2])