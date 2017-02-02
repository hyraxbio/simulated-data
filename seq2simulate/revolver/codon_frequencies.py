import numpy

class CodonTable(object):
    """
    A simple implementation of the standard genetic code table.

    Args:
        stop_codons: include stop codons (default True)
    """
    def __init__(self, stop_codons=True):
        bases = ['t', 'c', 'a', 'g']
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        self.codon_dict = dict(zip(codons, amino_acids))
        if not stop_codons:
            for codon in ['taa', 'tag', 'tga']:
                self.codon_dict.pop(codon)        
        for codon,aa in self.codon_dict.items():
            setattr(self, codon, aa)

    def __str__(self):
        return '<CodonTable {}>'.format(''.join(sorted([aa for aa in self.codon_dict.values()])))

def FEqual(*args, **kwargs):
    """
    Returns a dictionary of equal codon frequencies.

    Args:
        codon_table: an instance of CodonTable
    """
    try:
        codon_table = kwargs['codon_table'] 
    except:
        codon_table = None
    if codon_table is None:
        codon_table = CodonTable(stop_codons=False)
    codon_freq = {}
    for i in codon_table.codon_dict:
        codon_freq[i] = 1.0/len(codon_table.codon_dict)
    return codon_freq
    
def F1x4(sequence, codon_table=None):
    """
    Estimates equilibrium codon frequencies from sequence as in Goldman (1994).
    Codon frequency is the product of the frequencies of the three nucleotides in
    the triplet.

    Args:
        sequence: an instance of models.Sequence
        codon_table: an instance of CodonTable

    Refs:
        Goldman, N., & Yang, Z. (1994). A codon-based model of nucleotide
        substitution for protein-coding DNA sequences. Molecular Biology and Evolution,
        11(5), 725-736.
    """
    if codon_table is None:
        codon_table = CodonTable(stop_codons=False)
    seq = sequence.seq
    nfreq = {i:seq.count(i)/float(len(seq)) for i in 'atgc'}
    codons = sorted(codon_table.codon_dict.keys())
    codon_freq = []
    for codon in codons:
        codon_freq.append(numpy.product([nfreq[n] for n in codon]))
    codon_freq = numpy.array(codon_freq)
    codon_freq = codon_freq/codon_freq.sum()
    return dict(zip(codons, codon_freq))
