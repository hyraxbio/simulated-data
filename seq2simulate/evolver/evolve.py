import trees, models, codon_frequencies
import numpy

model_qfuncs = {
    'simple_goldman': models.goldman_Q,
}

def evolve_sequence(sequence, 
    t=0.1, 
    omega=1.0, 
    kappa=2.0, 
    codon_freq=None, 
    scale_q=True, 
    model='simple_goldman'):

    qfunc = model_qfuncs[model]
    q = qfunc(kappa=kappa, omega=omega, codon_freq=codon_freq, scale_q=scale_q, return_dict=False)
   
    loci = models.parse_sequence_to_loci(sequence)
    for locus in loci:
        for codon in locus.codons:
            codon.seq = models.call_mutation_from_q(codon.seq, q, t=t)
    return models.parse_loci_to_sequence_string(loci)

def evolve_sequence_with_q(sequence, q, t=0.1): 
    if not isinstance(t, (float, int)):
        raise ValueError('t must be a number')
    if not isinstance(q, numpy.ndarray):
        raise ValueError('q must be a NumPy array')
    loci = models.parse_sequence_to_loci(sequence)
    for locus in loci:
        for codon in locus.codons:
            codon.seq = models.call_mutation_from_q(codon.seq, q, t=t)
    return models.parse_loci_to_sequence_string(loci)

def evolve_tree(sequence,
    taxa=10,
    t=0.01, 
    omega=1.0, 
    kappa=2.0, 
    codon_freq=None, 
    scale_q=True, 
    model='simple_goldman'):

    """
    Evolve a parent DNA sequence into a set of daughter sequences (taxa) by:
        1. generating a random phylogenetic tree
        2. intantiating a mutational model (e.g. Goldman-Yang-like by default) represented by Q-matrix
        3. mutate sequence according to tree shape using model 

    Args:
        sequence: string of DNA nucleotides
        taxa: number of daughter sequences to evolve
        t: evolution time or branch length
        omega: dN/dS 
        kappa: ratio of transition to transversion rates
        codon_freq: dictionary of codon_frequencies, also know as equilibrium
        frequencies

        scale_q: scales Q so that the average rate of substitution at
        equilibrium equals 1. Branch lengths are thus expected number of nucleotide
        substitutions per codon. See Goldman (1994).
        model: mutational model, 'simple_goldman' will use a Goldman-Yang-like model

    Returns:
        tree instance populated with new sequence strings
    """
    sequence = sequence.lower()
    qfunc = model_qfuncs[model]
    q = qfunc(kappa=kappa, omega=omega, codon_freq=codon_freq, scale_q=scale_q, return_dict=False)
   
    tree = trees.random_tree(taxa)
    tree.value = sequence
    for node in trees.get_list_of_tree_nodes(tree)[1:]:
        node.value = evolve_sequence_with_q(node.parent.value, q, t=t)
    return tree 
        
def evolve(sequence,
    taxa=10,
    t=0.01, 
    omega=1.0, 
    kappa=2.0, 
    codon_freq=None, 
    scale_q=True, 
    model='simple_goldman'):

    tree = evolve_tree(**locals())
    return [i.value for i in trees.get_list_of_tree_leaves(tree)]


def print_mutations(old_sequence, new_sequence, colour=True):
    """
    Simple string comparison.
    """

    sim = ''
    for i,j in zip(old_sequence, new_sequence):
        if i==j:
            sim += ' '
        else:
            sim += '|'

    old_sequence_aa = [i for i in models.convert_dna_to_aa(old_sequence)] 
    new_sequence_aa = [i for i in models.convert_dna_to_aa(new_sequence)]
  
    if colour:
        red = ['\x1b[31m', '\x1b[0m']
        for i in range(len(old_sequence_aa)):
            if old_sequence_aa[i] != new_sequence_aa[i]:
                new_sequence_aa[i] = red[0]+new_sequence_aa[i]+red[1]
                
    old_sequence_aa = '..'.join(old_sequence_aa)    
    new_sequence_aa = '..'.join(new_sequence_aa)    
              

    print(old_sequence_aa)
    print(old_sequence)
    print(sim)
    print(new_sequence)
    print(new_sequence_aa)
    print('\n')

if __name__=='__main__':
    pass 
