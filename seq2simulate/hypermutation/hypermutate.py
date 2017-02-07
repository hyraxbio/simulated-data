import re
import mutation_probabilities


class Sequence(object):
    def __init__(self, sequence):
        """
        Convenience class to deal with DNA sequences.
        """
        if not all(nucleotide in 'ATGC' for nucleotide in set(sequence)):
            raise ValueError('Sequences may only contain A, T, G, and/or C.')

        self.sequence = sequence

    def get_motif_indices(self, motif):
        """
        Find indices of substring motif in sequence.
    
        Args:
            motif: a substring
    
        Returns: a list of left indices
        """
        indices = [m.start() for m in re.finditer(re.escape(motif), self.sequence)]
        return indices

    def index_all_motifs(self, motif_probabilities=None):
        """
        Returns a dictionary of all motif indices in the sequence.

        Args:
            motif_probabilities: instance of mutation_probabilities
        """
        if motif_probabilities is None:
            motif_probabilities = mutation_probabilities.KijakProbabilities()

        motif_dict = {}
        for motif in motif_probabilities.motifs:
            motif_dict[motif] = self.get_motif_indices(motif) 

        return motif_dict


if __name__ == '__main__':
    with open('../../data/split/Seq6_Sus', 'r') as f: 
        seq = f.read()
    seq = seq.replace('\n', '')[9:]
    km = mutation_probabilities.KijakProbabilities()     
    pass
