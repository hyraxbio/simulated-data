import unittest
import hypermutate
import mutation_probabilities
import numpy

class HypermutTester(unittest.TestCase):


    def setUp(self):
        self.old_sequence = 'atgcaacggcgattatacgtatcgtgcatcgatcatcgcatgcaacggcgattatacgtatcgtgcatcgatcatcgc'.upper()
        self.seq = hypermutate.Sequence(self.old_sequence)


    def test_sequence_init(self):
        self.assertEqual(self.old_sequence, self.seq.sequence)
        self.old_sequence += 'X'
        with self.assertRaises(ValueError):
            self.seq = hypermutate.Sequence(self.old_sequence)


    def test_get_motif_indices(self):
        indices = self.seq.get_motif_indices('ATG')
        self.assertTrue(indices, [0, 39])

    def test_index_all_motifs(self):
        indices = self.seq.index_all_motifs()     

    def test_motif_probabilities(self):
        motif_probabilities = mutation_probabilities.MotifProbabilities() 
        self.assertTrue(numpy.isclose(sum(motif_probabilities.motifs.values()), 1))

    def test_kijak_probabilities(self):
        motif_probabilities = mutation_probabilities.KijakProbabilities() 
        self.assertTrue(numpy.isclose(sum(motif_probabilities.motifs.values()), 1))
