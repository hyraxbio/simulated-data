import random
import numpy
import re
import mutation_probabilities
random.seed()

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
        else:
            self.motif_probabilities = motif_probabilities

    def _get_motif_indices(self, motif):
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
            motif_indices = self._get_motif_indices(motif) 
            if len(motif_indices) > 0:
                motif_dict[motif] = motif_indices

        self.sequence_motif_dict = motif_dict

   
    @property 
    def _num_motifs(self):
        """
        Returns number of substrate motifs currently present in sequence.
        """
        try:
            _num_motifs = sum([len(j) for i, j in self.sequence_motif_dict.iteritems()])
        except AttributeError:
            _num_motifs = 0
        return _num_motifs

    @property
    def _motif_probs(self):
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
        Attempt to perform n mutations of self.sequence or until no more substrate motifs exist.

        Args:
            n: int
        """

        if self.sequence_motif_dict is None:
            self.index_all_motifs()
        if n > self._num_motifs:
            n = self._num_motifs
            #raise ValueError('n must be less than or equal to _num_motifs ({})'.format(self._num_motifs))
            
        mutations = 0
        while mutations < n:
            self.sequence = self._mutate(self.sequence)
            self.index_all_motifs()
            mutations += 1

    def _mutate(self, sequence): 
        """
        Perform a single random mutation at the leftmost G residue in a motif
        with probabilities defined by Sequence.motif_probabilites.
        """
        if self.sequence_motif_dict is None:
            self.index_all_motifs()
        motif_index = [i > random.uniform(0, 1) for i in sorted(self._motif_probs_cum)].index(True)
        current_prob = sorted(self._motif_probs_cum)[motif_index]
        current_motif = self._motif_probs_cum[current_prob]
        current_index = random.choice(self.sequence_motif_dict[current_motif])
        mutation_index = current_index + current_motif.index('G')
        return sequence[:mutation_index] + 'A' + sequence[1 + mutation_index:]

    @property 
    def _motif_probs_cum(self):    
        return self._get_cumulative_p(self._motif_probs)

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
    
def print_mutations(old_sequence, new_sequence, colour=True):
    """
    Simple string comparison.
    """

    sim = ''
    for i, j in zip(old_sequence, new_sequence):
        if i == j:
            sim += ' '
        else:
            sim += '|'

    print(old_sequence)
    print(sim)
    print(new_sequence)
    print('\n')
        
if __name__ == '__main__':
    pass
