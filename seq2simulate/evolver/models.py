from numpy import array, ones, identity
import numpy
from codon_frequencies import CodonFrequencies




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
    Locus object containing codon(s) for grouping mutations by a reference codon.
    
    Args:
        codons: a list of Codon objects or nucleotide strings (typically triplets).
    """

    def __init__(self, codons=[], loc=0):
        self.loc = loc
        self._codons = []
        for codon in codons:
            self.add_codon(codon)             
   
    @property 
    def loc_aa(self):
        return self.loc+1

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

    @property
    def sequence(self):
        """
        Returns compiled sequence string.
        """
        return ''.join([codon.seq for codon in self.codons])
        
    @property    
    def location_sequence(self):
        """
        Returns codon location string and sequence string.
        """
        sequence = [codon.seq for codon in self.codons]
        return [self.loc_aa]*len(sequence), sequence

class Sequence(object):
    """
    Convenience class to collect a set of loci.
    """
    def __init__(self, seq):
        self.loci = parse_sequence_to_loci(seq) 

    @property
    def codons(self):
        return [j for i in [locus.codons for locus in self.loci] for j in i]

    @property
    def seq(self):
        return parse_loci_to_sequence_string(self.loci)
    
    @property
    def location_seq(self):
        return parse_loci_to_sequence(self.loci)
    
def parse_sequence_to_loci(sequence):
    nfloor = len(sequence)//3.0
    nfloat = len(sequence)/3.0
    if nfloat > nfloor:
        nfloor += 1
    nfloor = int(nfloor)
    codon_strings = [sequence[i*3:i*3+3] for i in range(nfloor)]
    loci = [Locus(codons=[j], loc=i) for i,j in enumerate(codon_strings)]
    return loci

def parse_loci_to_sequence(loci):
    locations = []
    sequences = []
    for locus in loci:
        location, sequence = locus.location_sequence
        locations.append(location)
        sequences.append(sequence)
    return [j for i in locations for j in i], [j for i in sequences for j in i]

def parse_loci_to_sequence_string(loci):
    sequences = []
    return ''.join([locus.sequence for locus in loci])

def goldman_Q(kappa=2.0, omega=1.0, codon_freq_dict=None):
    ct = CodonTable(stop_codons=False)
    codons = sorted(ct.codon_dict)
    n = len(codons)
    q = ones([n, n])
    idmatrix = identity(n)
    cf = CodonFrequencies().paml_fcodon

    for i in range(n):
        codon1 = codons[i]
        for j in range(n):
            codon2 = codons[j] 
            q[i, j] = mutation_rate(codon1, codon2, codon_table=ct, codon_freq=cf)

    # diagonal equals negative sum of row
    q[idmatrix==1] = 0.0
    q[idmatrix==1] = -q.sum(1)

    # convert matrix to codon-indexed dict
    qdict = dict(zip(codons, [dict(zip(codons, qq)) for qq in q]))

    return q, qdict

def mutation_category(codon1, codon2, codon_table=None):
    """
    Args:
        codon1/codon2: three-letter strings

    Returns category of mutation one of:
        ['multisite'],
        ['synonymous', 'transversion],
        ['synonymous', 'transition],
        ['nonsynonymous', 'transversion],
        ['nonsynonymous', 'transition],
    """
    if not isinstance(codon1, str):
        raise ValueError('codon must be a three-letter str')
    if len(codon1) != 3:
        raise ValueError('codon must be a three-letter str')
    if not isinstance(codon2, str):
        raise ValueError('codon must be a three-letter str')
    if len(codon2) != 3:
        raise ValueError('codon must be a three-letter str')
    if not isinstance(codon_table, CodonTable):
        raise ValueError('codon must be a three-letter str')

    mscores = mutations(codon1, codon2)
    
    cat = None    

    if sum([nt > -1 for nt in mscores]) > 1:
        cat = ['multisite']
    else:
        nscore = max(mscores)

        nucleotide_mutation_type = None
        if nscore == 0:
            nucleotide_mutation_type = 'transition'
        if nscore == 1:
            nucleotide_mutation_type = 'transversion'

        aminoacid_mutation_type = None
        if getattr(codon_table, codon1) == getattr(codon_table, codon2):
            aminoacid_mutation_type = 'synonymous'
        else:
            aminoacid_mutation_type = 'nonsynonymous'

        cat = [aminoacid_mutation_type, nucleotide_mutation_type]

    return cat

def mutation_rate(codon1, codon2, 
    kappa=2.0, 
    omega=1.0, 
    codon_table=None, 
    codon_freq=None
    ):
    """
    Args:
        codon1/codon2: three-letter strings
        codon_table: CodonTable instance
        kappa: transition/transversion ratio
        omega: dN/dS
        codon_freq: dictionary of codon equilibrium frequencies

    Returns instantaneous rate of mutation using the simplified Goldman model in Nielsen and Yang (1998).
    """
    if not isinstance(codon1, str):
        raise ValueError('codon must be a three-letter str')
    if len(codon1) != 3:
        raise ValueError('codon must be a three-letter str')
    if not isinstance(codon2, str):
        raise ValueError('codon must be a three-letter str')
    if len(codon2) != 3:
        raise ValueError('codon must be a three-letter str')
    if not isinstance(codon_table, CodonTable):
        raise ValueError('codon must be a three-letter str')
    if not isinstance(omega, (int, float)):
        raise ValueError('omega must be a number')
    if not isinstance(kappa, (int, float)):
        raise ValueError('kappa must be a number')
    if not isinstance(codon_freq, dict):
        raise ValueError('codon_freq must be a dictionary')


    cat = mutation_category(codon1, codon2, codon_table=codon_table)
    rate = 0
    if cat == 'multisite':
        rate = 0 
    if cat == ['synonymous', 'transversion' ]:
        rate = codon_freq[codon2]
    if cat == ['synonymous', 'transition' ]:
        rate = codon_freq[codon2]*kappa
    if cat == ['nonsynonymous', 'transversion']:
        rate = codon_freq[codon2]*omega
    if cat == ['nonsynonymous', 'transition']:
        rate = codon_freq[codon2]*kappa*omega
    return rate 

def diff_index(s1, s2):
    index = []
    for i,j in zip(s1, s2):
        if i==j:
            index.append(0)
        else: 
            index.append(1)
    return index
     
def mutations(seq1, seq2):
    """
    Returns list of mutation values: 0 for transversion, 1 for transition, -1 for neither.

    Args:
        seq1, seq2: nucleotide strings
    """
    mut_dict = {
        'at': 1, 
        'ag': 0,
        'ac': 1,
        'ta': 1,
        'tg': 1,
        'tc': 0,
        'ga': 0,
        'gt': 1,
        'gc': 1,
        'ca': 1,
        'ct': 0,
        'cg': 1,
    }
    mutations = []
    for i,j in zip(seq1, seq2):
        mutations.append(mut_dict.get(i+j, -1))
    return mutations
            

if __name__=='__main__':
    q, qdict = goldman_Q()
    pass

