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
        new_sequence = evolve.evolve_sequence_with_q(self.old_sequence, q, t=1)
        self.assertNotEqual(self.old_sequence, new_sequence) # this has a very low probability of failing

    def test_evolve_tree(self):
        t10 = evolve.evolve_tree(self.old_sequence, taxa=10, t=0.1, omega=1.1, kappa=1.5)
        new_sequences = [i.value for i in trees.get_list_of_tree_leaves(t10)]
        self.assertEqual(len(new_sequences), 10)
        for i in new_sequences:
            self.assertNotEqual(self.old_sequence, i)

if __name__=='__main__':
    unittest.main()
