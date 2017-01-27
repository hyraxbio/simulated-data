
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

FEqual = {
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

