import unittest
import hypermutate
import mutation_probabilities
import numpy

class HypermutTester(unittest.TestCase):


    def setUp(self):
        self.old_sequence = 'atgcaacggcgattatacgtatcgtgcatcgatcatcgcatgcaacggcgattatacgtatcgtgcatcgatcatggtcgc'.upper()
        self.seq = hypermutate.Sequence(self.old_sequence)


    def test_sequence_init(self):
        self.assertEqual(self.old_sequence, self.seq.sequence)
        self.old_sequence += 'X'
        with self.assertRaises(ValueError):
            self.seq = hypermutate.Sequence(self.old_sequence)


    def test_get_motif_indices(self):
        indices = self.seq.get_motif_indices('ATG')
        self.assertEqual(indices, [0, 39, 73])


    def test_index_all_motifs(self):
        self.seq.index_all_motifs()
        self.assertEqual(self.seq.sequence_motif_dict['CGGC'], [6, 45])


    def test_motif_probabilities(self):
        motif_probabilities = mutation_probabilities.MotifProbabilities() 
        self.assertTrue(numpy.isclose(sum(motif_probabilities.motifs.values()), 1))


    def test_kijak_probabilities(self):
        motif_probabilities = mutation_probabilities.KijakProbabilities() 
        self.assertTrue(numpy.isclose(sum(motif_probabilities.motifs.values()), 1))


    def test_num_motifs(self):
        self.seq.index_all_motifs()
        self.assertEqual(self.seq.num_motifs, 7)


    def test_motif_probs(self):
        self.seq.index_all_motifs()
        expected_answer = {'CGGC': 0.0, 'TGGT': 0.3890518084067808, 'GA': 0.6109481915932192}
        for i, j in expected_answer.iteritems():
            self.assertTrue(numpy.isclose(self.seq.motif_probs[i], j))
