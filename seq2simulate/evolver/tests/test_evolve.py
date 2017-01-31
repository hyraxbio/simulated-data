import unittest
import models
import codon_frequencies
import trees
import evolve
from numpy import isclose, identity

class EvolveTester(unittest.TestCase):

    def setUp(self):
        self.old_sequence = 'atgcaacggcgattatacgtatcgtgcatcgatcatcgcatgcaacggcgattatacgtatcgtgcatcgatcatcgc'

    def test_evolve_sequence_with_q(self):
        qfunc = evolve.model_qfuncs['simple_goldman']
        q = qfunc(scale_q=True)
        new_sequence0 = evolve.evolve_sequence_with_q(models.Sequence(self.old_sequence), q, t=1)
        new_sequence1 = evolve.evolve_sequence_with_q(models.Sequence(self.old_sequence), q, t=1)
        new_seq0 = new_sequence0.seq
        new_seq1 = new_sequence1.seq

        self.assertNotEqual(self.old_sequence, new_seq0) # this has a very low probability of failing
        self.assertNotEqual(self.old_sequence, new_seq1) # this has a very low probability of failing
        self.assertNotEqual(new_seq0, new_seq1)

    def test_evolve_tree(self):
        t10 = evolve.evolve_tree(models.Sequence(self.old_sequence), taxa=10, t=0.1, omega=1.1, kappa=1.5)
        new_sequences = [i.value for i in trees.get_list_of_tree_leaves(t10)]
        self.assertEqual(len(new_sequences), 10)
        for i in new_sequences:
            self.assertNotEqual(self.old_sequence, i.seq) # this has a very low probability of failing

    def test_evolve(self):
        new_sequences = evolve.evolve(self.old_sequence, taxa=10, t=0.1, omega=1.1, kappa=1.5, lmbda=0.01, ti_td=0.1)
        self.assertEqual(len(new_sequences), 10)
        for i in new_sequences:
            self.assertNotEqual(self.old_sequence, i) # this has a very low probability of failing

    def test_evolve_no_insertion(self):
        new_sequences = evolve.evolve(self.old_sequence, taxa=10, t=0.1, omega=1.1, kappa=1.5, lmbda=0.01, ti_td=0.0)
        self.assertEqual(len(new_sequences), 10)
        for i in new_sequences:
            self.assertEqual(len(self.old_sequence), len(i))

    def test_evolve_with_high_insertion(self):
        new_sequences = evolve.evolve(self.old_sequence, taxa=10, t=0.1, omega=1.1, kappa=1.5, lmbda=0.1, ti_td=10.0)
        self.assertEqual(len(new_sequences), 10)
        len_new_sequences = sum([len(i) for i in new_sequences])
        self.assertTrue(len_new_sequences != 10*len(self.old_sequence))

    def test_evolve_returns_different_evolved_sequences(self):
        new_sequences0 = evolve.evolve(self.old_sequence, taxa=10, t=0.1, omega=1.1, kappa=1.5, lmbda=0.01, ti_td=0.1)
        new_sequences1 = evolve.evolve(self.old_sequence, taxa=10, t=0.1, omega=1.1, kappa=1.5, lmbda=0.01, ti_td=0.1)
        for i,j in zip(new_sequences0, new_sequences1):
            self.assertNotEqual(i, j)

    def test_compile_histories(self):
        t10 = evolve.evolve_tree(models.Sequence(self.old_sequence), taxa=10, t=0.1, omega=1.1, kappa=1.5)
        histories = evolve.compile_histories(t10)
        self.assertIsInstance(histories, dict)
        for i in trees.get_list_of_tree_nodes(t10):
            self.assertIn(i.id, histories)

    def test_evolve_with_log(self):
        new_sequences, mutations = evolve.evolve(self.old_sequence, taxa=10, t=0.1, omega=1.1, kappa=1.5, lmbda=0.01, ti_td=0.1, log=True)
        print(mutations)

if __name__=='__main__':
    unittest.main()
