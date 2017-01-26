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
        self.assertEqual(l.loc, 0)
        self.assertEqual(l.loc_aa, 1)
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
        
    def test_parse_hanging_sequence_to_loci(self):
        s = 'atgatgccagtcgatcgatcgtagcatcgtagctgtagca'
        l = models.parse_sequence_to_loci(s)
        for i,seq in enumerate(['atg', 'atg', 'cca', 'gtc', 'gat', 'cga', 'tcg', 'tag', 'cat', 'cgt', 'agc', 'tgt', 'agc','a--']):
            self.assertEqual(len(l[i].codons), 1)
            self.assertEqual(l[i].codons[0].seq, seq)

    def test_parse_sequence_to_loci(self):
        s = 'atgatgccagtcgatcgatcgtagcatcgtagctgtagc'
        l = models.parse_sequence_to_loci(s)
        for i,seq in enumerate(['atg', 'atg', 'cca', 'gtc', 'gat', 'cga', 'tcg', 'tag', 'cat', 'cgt', 'agc', 'tgt', 'agc']):
            self.assertEqual(len(l[i].codons), 1)
            self.assertEqual(l[i].codons[0].seq, seq)

    def test_parse_loci_to_sequence(self):
        s0 = 'atgatgccagtcgatcgatcgtagcatcgtagctgtagca'
        l = models.parse_sequence_to_loci(s0)
        l[1].add_codon('ttt')
        l1, s1 = models.parse_loci_to_sequence(l)
        self.assertEqual(l1, [1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
        self.assertEqual(s1, ['atg', 'atg', 'ttt', 'cca', 'gtc', 'gat', 'cga', 'tcg', 'tag', 'cat', 'cgt', 'agc', 'tgt', 'agc', 'a--'])

    def test_parse_loci_to_sequence_string(self):
        s0 = 'atgatgccagtcgatcgatcgtagcatcgtagctgtagc'
        l = models.parse_sequence_to_loci(s0)
        s1 = models.parse_loci_to_sequence_string(l)
        self.assertEqual(s1, s0)


if __name__=='__main__':
    unittest.main()
