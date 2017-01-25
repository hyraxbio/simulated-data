import unittest
import models

class ModelTester(unittest.TestCase):
  
    def test_codontable_init(self):
        ct = models.CodonTable()
        self.assertEqual(ct.__str__(), '<CodonTable ***AAAACCDDEEFFGGGGHHIIIKKLLLLLLMNNPPPPQQRRRRRRSSSSSSTTTTVVVVWYY>')

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
        self.assertEqual(c.__str__(), '<Codon --->')

    def test_codon_seq_validation(self):
        c = models.Codon()
        for seq in [1, 'abc']:
            with self.assertRaises(ValueError):
                c.seq = seq
        c.seq = 'a'
        self.assertEqual(c.seq, 'a--')

    def test_locus_init(self):
        l = models.Locus()
        self.assertEqual(l.__str__(), '<Locus >')
        l = models.Locus(codons=[models.Codon(), 'atg'])
        self.assertEqual(len(l.codons), 2)
        self.assertEqual(l.codons[0].seq, '---')
        self.assertEqual(l.codons[1].seq, 'atg')

    def test_locus_add_codon(self):
        c = models.Codon()
        l = models.Locus()
        for codon in [c, 'atg']:
            l.add_codon(codon)
        self.assertEqual(len(l.codons), 2)
        self.assertEqual(l.codons[0].seq, '---')
        self.assertEqual(l.codons[1].seq, 'atg')

    def test_locus_cannot_add_spurious_codons(self):
        c = models.Codon()
        with self.assertRaises(ValueError):
            l = models.Locus(codons=[c, 'abc'])
        

if __name__=='__main__':
    unittest.main()
