import numpy
import re
import mutation_probabilities


class Sequence(object):

    def __init__(self, sequence, motif_probabilities=None):
        """
        Convenience class to deal with DNA sequences.

        Args:
            sequence: a DNA string
            motif_probabilities: instance of mutation_probabilities
        """
        if not all(nucleotide in 'ATGC' for nucleotide in set(sequence)):
            raise ValueError('Sequences may only contain A, T, G, and/or C.')

        self.sequence = sequence
        self.sequence_motif_dict = None
        if motif_probabilities is None:
            self.motif_probabilities = mutation_probabilities.KijakProbabilities()

    def get_motif_indices(self, motif):
        """
        Find indices of substring motif in sequence.
    
        Args:
            motif: a substring
    
        Returns: a list of left indices
        """
        indices = [m.start() for m in re.finditer(re.escape(motif), self.sequence)]
        return indices

    def index_all_motifs(self):
        """
        Returns a dictionary of all motif indices in the sequence.
        """

        motif_dict = {}
        for motif in self.motif_probabilities.motifs:
            motif_indices = self.get_motif_indices(motif) 
            if len(motif_indices) > 0:
                motif_dict[motif] = motif_indices

        self.sequence_motif_dict = motif_dict

   
    @property 
    def num_motifs(self):
        try:
            num_motifs = sum([len(j) for i, j in self.sequence_motif_dict.iteritems()])
        except AttributeError:
            num_motifs = 0
        return num_motifs

    @property
    def motif_probs(self):
        """
        A shortened rescaled dictionary of motif mutation probabilities
        including only the present motifs.

        """
        try:
            remaining_probabilities_sum = sum([self.motif_probabilities.motifs[i] for i in self.sequence_motif_dict])
            remaining_probabilities = {i: self.motif_probabilities.motifs[i] / remaining_probabilities_sum for i in self.sequence_motif_dict}
        except AttributeError:
            remaining_probabilities = {}
        return remaining_probabilities
            

    def mutate_sequence(self, n):
        """
        Perform n mutations of self.sequence.

        Args:
            n: int
        """
        if self.sequence_motif_dict is None:
            print('First indexing sequence.')
            self.index_all_motifs()

        mutations = []
        #while len(mutations) < n:
   
    @property 
    def motif_probs_cum(self):    
        return self._get_cumulative_p(self.motif_probs)

    def _get_cumulative_p(self, p):
        """
        Returns cumulatively probabilities from probability dict sorted by
        cumulative probability. This function is useful for sampling mutations.
    
        Args:
            p: sorted cumulative probability dict
        """
     
        motifs = p.keys()
        probabilities = p.values()
        p_cumsum = numpy.array(probabilities).cumsum()
        return dict(zip(p_cumsum, motifs))
    
        
if __name__ == '__main__':
    seq = 'atgcaacggcgattatacgtatcgtgcatcgatcatcgcatgcaacggcgattatacgtatcgtgcatcgatcatggtcgc'.upper()
    seq = Sequence(seq)
    #km = mutation_probabilities.KijakProbabilities()     
    i = seq.index_all_motifs()

    pass
