import unittest
import models

class ModelTester(unittest.TestCase):
  
    def test_codontable_init(self):
        ct = models.CodonTable()
        self.assertEqual(ct.__repr__(), '<CodonTable ***AAAACCDDEEFFGGGGHHIIIKKLLLLLLMNNPPPPQQRRRRRRSSSSSSTTTTVVVVWYY>')

    def test_codontable_can_exclude_stop_codons(self):
        cts = models.CodonTable(stop_codons=False)
        for codon in ['taa', 'tag', 'tga']:
            self.assertNotIn(codon, cts.codon_dict)
            with self.assertRaises(AttributeError):
                getattr(cts, codon)
            

if __name__=='__main__':
    unittest.main()
