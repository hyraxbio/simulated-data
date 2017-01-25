import unittest
import models

class ModelTester(unittest.TestCase):
  
    def test_codontable_init(self):
        ct = models.CodonTable()
        self.assertEqual(ct.__repr__(), '<CodonTable ***AAAACCDDEEFFGGGGHHIIIKKLLLLLLMNNPPPPQQRRRRRRSSSSSSTTTTVVVVWYY>')

    def test_codontable_can_exclude_stop_codons(self):
        ct = models.CodonTable()
        cts = models.CodonTable(stop_codons=False)
        for codon in ['taa', 'tag', 'tga']:
            self.assertIn(codon, ct.codon_dict)
            getattr(ct, codon)
            self.assertNotIn(codon, cts.codon_dict)
            with self.assertRaises(AttributeError):
                getattr(cts, codon)
            
    def test_codon_init(self):
        c = models.Codon()
        self.assertEqual(c.__repr__(), '<Codon --->')

    def test_codon_seq_validation(self):
        c = models.Codon()
        for seq in [1, 'abc']:
            with self.assertRaises(ValueError):
                c.seq = seq
        c.seq = 'a'
        self.assertEqual(c.seq, 'a--')

if __name__=='__main__':
    unittest.main()
