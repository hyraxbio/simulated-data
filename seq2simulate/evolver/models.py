from numpy import array, ones, identity
import numpy





class ValidationMixin(object):

    def _is_valid_dna(self, seq):
        valids = 'ATGCatgc-'
        for i in set(seq):
            if i not in valids: 
                return False
        return True, ''


codon_frequencies = {
    'aaa': 0.016393442623,
    'aac': 0.016393442623,
    'aag': 0.016393442623,
    'aat': 0.016393442623,
    'aca': 0.016393442623,
    'acc': 0.016393442623,
    'acg': 0.016393442623,
    'act': 0.016393442623,
    'aga': 0.016393442623,
    'agc': 0.016393442623,
    'agg': 0.016393442623,
    'agt': 0.016393442623,
    'ata': 0.016393442623,
    'atc': 0.016393442623,
    'atg': 0.016393442623,
    'att': 0.016393442623,
    'caa': 0.016393442623,
    'cac': 0.016393442623,
    'cag': 0.016393442623,
    'cat': 0.016393442623,
    'cca': 0.016393442623,
    'ccc': 0.016393442623,
    'ccg': 0.016393442623,
    'cct': 0.016393442623,
    'cga': 0.016393442623,
    'cgc': 0.016393442623,
    'cgg': 0.016393442623,
    'cgt': 0.016393442623,
    'cta': 0.016393442623,
    'ctc': 0.016393442623,
    'ctg': 0.016393442623,
    'ctt': 0.016393442623,
    'gaa': 0.016393442623,
    'gac': 0.016393442623,
    'gag': 0.016393442623,
    'gat': 0.016393442623,
    'gca': 0.016393442623,
    'gcc': 0.016393442623,
    'gcg': 0.016393442623,
    'gct': 0.016393442623,
    'gga': 0.016393442623,
    'ggc': 0.016393442623,
    'ggg': 0.016393442623,
    'ggt': 0.016393442623,
    'gta': 0.016393442623,
    'gtc': 0.016393442623,
    'gtg': 0.016393442623,
    'gtt': 0.016393442623,
    'tac': 0.016393442623,
    'tat': 0.016393442623,
    'tca': 0.016393442623,
    'tcc': 0.016393442623,
    'tcg': 0.016393442623,
    'tct': 0.016393442623,
    'tgc': 0.016393442623,
    'tgg': 0.016393442623,
    'tgt': 0.016393442623,
    'tta': 0.016393442623,
    'ttc': 0.016393442623,
    'ttg': 0.016393442623,
    'ttt': 0.016393442623,
}



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

def goldman_Q():
    ct = CodonTable(stop_codons=False)
    codons = sorted(ct.codon_dict)
    n = len(codons)
    q = ones([n, n])
    idmatrix = identity(n)

    pars = {
        'k': 0.2,
        'w': 1.1,
    }


    #for i in range(n):
    #    codon1 = codons[i]
    #    for j in range(n):
    #        codon2 = codons[j] 

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

        cat = [nucleotide_mutation_type, aminoacid_mutation_type]

    return cat

def mutation_rate(codon1, codon2, codon_table=None):
    """
    Args:
        codon1/codon2: three-letter strings

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

    cat = mutation_category(codon1, codon2, codon_table=codon_table)
    rate = 0
    if cat == 'multisite':
        rate = 0 
    if cat == ['transversion', 'synonymous']:
        rate = codon_frequencies[codon2]
    if cat == ['transition', 'synonymous']:
        rate = codon_frequencies[codon2]
    if cat == ['transversion', 'nonsynonymous']:
        rate = codon_frequencies[codon2]
    if cat == ['transition', 'nonsynonymous']:
        rate = codon_frequencies[codon2]
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
    #q, qdict = goldman_Q()
    codon_table = CodonTable(stop_codons=False)
    pass

