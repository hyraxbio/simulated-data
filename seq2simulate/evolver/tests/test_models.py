import unittest
import models
import codon_frequencies
from numpy import isclose, identity

class ModelTester(unittest.TestCase):
  
    def test_codontable_init(self):
        ct = codon_frequencies.CodonTable()
        self.assertEqual(ct.__str__(), '<CodonTable ***AAAACCDDEEFFGGGGHHIIIKKLLLLLLMNNPPPPQQRRRRRRSSSSSSTTTTVVVVWYY>')

    def test_codontable_can_exclude_stop_codons(self):
        ct = codon_frequencies.CodonTable()
        cts = codon_frequencies.CodonTable(stop_codons=False)
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

    def test_codon_delete(self):
        c = models.Codon(seq='atg')
        self.assertEqual(c.seq, 'atg')
        c.delete()
        self.assertEqual(c.seq, '---')

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

    def test_mutation_category_validation(self):
        ct = codon_frequencies.CodonTable(stop_codons=False)
        with self.assertRaises(ValueError):
            models.mutation_category('atg', 'acc', codon_table=[123])
        with self.assertRaises(ValueError):
            models.mutation_category('atg', 2, codon_table=ct)
        
    def test_mutation_category(self):
        ct = codon_frequencies.CodonTable(stop_codons=False)
        self.assertEqual(models.mutation_category('atg', 'acc', codon_table=ct), ['multisite'])
        self.assertEqual(models.mutation_category('gca', 'gcg', codon_table=ct), ['synonymous',    'transition'])
        self.assertEqual(models.mutation_category('gca', 'gct', codon_table=ct), ['synonymous',    'transversion'])
        self.assertEqual(models.mutation_category('atg', 'ata', codon_table=ct), ['nonsynonymous', 'transition'])
        self.assertEqual(models.mutation_category('atg', 'atc', codon_table=ct), ['nonsynonymous', 'transversion'])

    def test_mutation_rate_validation(self):
        ct = codon_frequencies.CodonTable(stop_codons=False)
        cf = codon_frequencies.FEqual
        with self.assertRaises(ValueError):
            models.mutation_rate('atg', 'acc')
        with self.assertRaises(ValueError):
            models.mutation_rate('atg', 'acc', codon_table=ct)
        with self.assertRaises(ValueError):
            models.mutation_rate('atg', 'acc', codon_freq=cf)
        with self.assertRaises(ValueError):
            models.mutation_rate('atg', 2, codon_table=ct, codon_freq=cf)
        
    def test_mutation_rate(self):
        ct = codon_frequencies.CodonTable(stop_codons=False)
        cf = codon_frequencies.FEqual
        self.assertEqual(models.mutation_rate('atg', 'acc', codon_table=ct, codon_freq=cf), 0)
        self.assertIsInstance(models.mutation_rate('gca', 'gcg', codon_table=ct, codon_freq=cf), float)

    def test_goldman_Q_validation(self):
        cf = codon_frequencies.FEqual
        with self.assertRaises(ValueError):
            models.goldman_Q(kappa='2')
        with self.assertRaises(ValueError):
            models.goldman_Q(omega='2')
        with self.assertRaises(ValueError):
            models.goldman_Q(codon_freq=['test'])
        
    def test_goldman_Q(self):
        cf = codon_frequencies.FEqual
        q, qdict = models.goldman_Q(codon_freq=cf, scale_q=False, return_dict=True)
        self.assertEqual(q.shape, (61, 61))
        for row in q:
            self.assertTrue(isclose(row.sum(), 0))

    def test_convert_q_to_p_validation(self):
        q = models.goldman_Q()
        with self.assertRaises(ValueError):
            models.convert_q_to_p(['blah'])
        with self.assertRaises(ValueError):
            models.convert_q_to_p(q, t='3')

    def test_convert_q_to_p(self):
        q = models.goldman_Q()
        p = models.convert_q_to_p(q, t=0)
        self.assertEqual(p.shape, (61, 61))
        self.assertTrue((p==identity(61)).all())
        p = models.convert_q_to_p(q, t=1)
        self.assertTrue((p.max() <= 1.0))

    def test_get_cumulative_p_validation(self):
        q = models.goldman_Q(scale_q=False)
        p = models.convert_q_to_p(q, t=10)
        with self.assertRaises(ValueError):
            models.convert_q_to_p(['blah'])
        with self.assertRaises(ValueError):
            models.convert_q_to_p(p, codon_table='asdf')

    def test_get_cumulative_p(self):
        q = models.goldman_Q(scale_q=False)
        p = models.convert_q_to_p(q, t=10)
        pc, pcod = models.get_cumulative_p(p)
        for i in pc:
            self.assertTrue(isclose(i[-1], 1))
        
    def test_get_mutation_from_cumulative_p_validation(self):
        q = models.goldman_Q(scale_q=False)
        p = models.convert_q_to_p(q, t=10)
        pc, pcod, pcdict = models.get_cumulative_p(p, return_dict=True)
        with self.assertRaises(ValueError):
            models.get_mutation_from_cumulative_p(3, pcdict)
        with self.assertRaises(ValueError):
            models.get_mutation_from_cumulative_p('aaa', 'blah')

    def test_get_mutation_from_cumulative_p(self):
        q = models.goldman_Q(scale_q=False)
        p = models.convert_q_to_p(q, t=10)
        pc, pcod, pcdict = models.get_cumulative_p(p, return_dict=True)
        old_codon = 'aaa'
        new_codon = models.get_mutation_from_cumulative_p(old_codon, pcdict)
                
    def test_call_sub_from_q_validation(self):
        q = models.goldman_Q(scale_q=False)
        with self.assertRaises(ValueError):
            models.call_sub_from_q(123, q, t=10)
        with self.assertRaises(ValueError):
            models.call_sub_from_q('aaa', 'blha', t=10)
        with self.assertRaises(ValueError):
            models.call_sub_from_q('aaa', q, t='s')

    def test_call_sub_from_q(self):
        q = models.goldman_Q(scale_q=False)
        old_codon = 'aaa'
        new_codon = models.call_sub_from_q(old_codon, q, t=10)
        self.assertIsInstance(new_codon, str)
        self.assertEqual(len(new_codon), 3)
        allowed_letters = 'atgc'
        for i in new_codon:
            self.assertIn(i, allowed_letters)

    def test_sample_model_mutation_probabilities_validation(self):
        q = models.goldman_Q(scale_q=False)
        with self.assertRaises(ValueError):
            models.sample_model_mutation_probabilities(123, q)
        with self.assertRaises(ValueError):
            models.sample_model_mutation_probabilities('ata', 'blah')
        with self.assertRaises(ValueError):
            models.sample_model_mutation_probabilities('ata', q, t=-1)
        with self.assertRaises(ValueError):
            models.sample_model_mutation_probabilities('ata', q, n=0)

    def test_sample_model_mutation_probabilities_validation(self):
        q = models.goldman_Q(scale_q=False)
        sample = models.sample_model_mutation_probabilities('aaa', q)
        self.assertIsInstance(sample, list)
        self.assertEqual(len(sample), 100) 

    def test_make_indel_force_deletion(self):
        c = models.Codon(seq='atg')
        l = models.Locus(codons=[c])
        models.make_indel(l, index=0, ti_td=0.0, codon_freq=None)
        self.assertIn('---', [i.seq for i in l.codons])
        
    def test_make_indel_force_insertion(self):
        c = models.Codon(seq='atg')
        l = models.Locus(codons=[c])
        models.make_indel(l, index=0, ti_td=1e12, codon_freq=None)
        self.assertEqual(len(l.codons), 2)
        self.assertFalse(any(i.seq == '---' for i in l.codons))
        



if __name__=='__main__':
    unittest.main()
