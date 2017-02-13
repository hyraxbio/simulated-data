import unittest
import hypermutate
import mutation_probabilities
import numpy

class HypermutTester(unittest.TestCase):


    def setUp(self):
        self.old_sequence = 'atgcaacggcgattatacgtatcgtgcatcgatcatcgcatgcaacggcgattatacgtatcgtgcatcgatcatggtcgc--'.upper()
        self.seq = hypermutate.Sequence(self.old_sequence)
    

    def test_sequence_init(self):
        self.assertEqual(self.old_sequence, self.seq.sequence)
        self.old_sequence += 'X'
        with self.assertRaises(ValueError):
            self.seq = hypermutate.Sequence(self.old_sequence)


    def test_get_motif_indices(self):
        indices = self.seq._get_motif_indices('ATG')
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


    def test__num_motifs(self):
        self.seq.index_all_motifs()
        self.assertEqual(self.seq._num_motifs, 7)


    def test__motif_probs(self):
        self.seq.index_all_motifs()
        expected_motifs = ['CGGC', 'TGGT', 'GA']
        expected_probs = [1e-18, 0.3890518084067808, 0.6109481915932192]
        for i, j in zip(expected_motifs, expected_probs):
            self.assertTrue(numpy.isclose(self.seq._motif_probs[i], j))


    def test__mutate(self):
        oldseq = self.seq.sequence
        self.seq._mutate()
        newseq = self.seq.sequence
        self.assertNotEqual(oldseq, newseq)
        diffs = [i for i in xrange(len(oldseq)) if oldseq[i] != newseq[i]]
        self.assertEqual(len(diffs), 1)
        self.assertEqual(oldseq[diffs[0]], 'G')
        self.assertEqual(newseq[diffs[0]], 'A')


    def test_mutate_sequence(self):
        oldseq = self.seq.sequence
        self.seq.mutate_sequence(3)
        newseq = self.seq.sequence
        diffs = [i for i in xrange(len(oldseq)) if oldseq[i] != newseq[i]]
        self.assertEqual(len(diffs), 3)
       
 
    def test_mutate_sequence_is_limited_to__num_motifs(self):
        oldseq = self.seq.sequence
        n = self.seq._num_motifs 
        self.seq.mutate_sequence(2*n)
        newseq = self.seq.sequence
        diffs = [i for i in xrange(len(oldseq)) if oldseq[i] != newseq[i]]
        self.assertEqual(len(diffs), n)


    def test_mutate_sequences(self):
        oldseqs = [self.seq.sequence]*5
        newseqs = hypermutate.mutate_sequences(oldseqs, 3)
        for i,j in zip(oldseqs, newseqs):
            diffs = [x for x in xrange(len(i)) if i[x] != j[x]]
            self.assertEqual(len(diffs), 3)


    def test_mutate_sequences_list_of_ns(self):
        oldseqs = [self.seq.sequence]*5
        newseqs = hypermutate.mutate_sequences(oldseqs, [3]*5)
        for i,j in zip(oldseqs, newseqs):
            diffs = [x for x in xrange(len(i)) if i[x] != j[x]]
            self.assertEqual(len(diffs), 3)



