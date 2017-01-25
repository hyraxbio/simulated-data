
class ValidationMixin(object):

    def _is_valid_dna(self, seq):
        valids = 'ATGCatgc-'
        for i in set(seq):
            if i not in valids: 
                return False
        return True, ''

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

class Codon(ValidationMixin, object):
    """
    Codon object containing triplet nucleotides (and possibly indels).
    
    Args:
        seq: the nucleotide sequence of the codon (default '---')
    """
    def __init__(self, seq='---'): 
        self.seq = seq

    @property
    def seq(self):
        return self._seq

    @seq.setter
    def seq(self, seq):
        if not isinstance(seq, str):
            raise ValueError('seq must be a str')
        #if len(seq) > 3:
        #    seq = seq[:3]    
        while len(seq) < 3:
            seq += '-'
        if not self._is_valid_dna(seq):
            raise ValueError('sequence is not valid')
        self._seq = seq
            
    def __str__(self):
        return '<Codon {}>'.format(self._seq)
        
class Locus(object):
    """
    Locus object containing codon(s) for indexing mutations by codon.
    
    Args:
        codons: a list of Codon objects or nucleotides
    """

    def __init__(self, codons=[]):
        self._codons = []
        for codon in codons:
            self.add_codon(codon)             
   
    @property 
    def codons(self):
        return self._codons

    @codons.setter 
    def codons(self, codons):
        if not isinstance(codons, list):
            raise ValueError('codons must be a list') 
        for codon in codons:
            if not isinstance(codon, Codon):
                raise ValueError('codons must contain Codons')
        self._codons = codons

    def add_codon(self, codon):
        if isinstance(codon, Codon):
            self._codons.append(codon)
        elif isinstance(codon, str):
            self._codons.append(Codon(seq=codon))
            
    def del_codon(self):
        self._codons = self._codons[:-1]

    def __str__(self):
        return '<Locus {}>'.format(' '.join([codon.__str__() for codon in self.codons]))

if __name__=='__main__':
    pass
