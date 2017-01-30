import trees, models, codon_frequencies
import numpy
from random import uniform, seed
from copy import deepcopy
seed()

model_qfuncs = {
    'simple_goldman': models.goldman_Q,
}

def evolve_sequence_with_q(sequence, q, t=0.1, lmbda=0.001, ti_td=0.1, indel_codon_freq=None):
    if not isinstance(sequence, models.Sequence):
        raise ValueError('sequence must be an instance of models.Sequence')
    if not isinstance(t, (float, int)):
        raise ValueError('t must be a number')
    if not isinstance(q, numpy.ndarray):
        raise ValueError('q must be a NumPy array')
    if not isinstance(lmbda, (float, int)):
        raise ValueError('lmbda must be a number')
    if not isinstance(ti_td, (float, int)):
        raise ValueError('ti_td must be a number')
    if not lmbda >= 0:
        raise ValueError('lmbda must be in range 0-1')
    if not lmbda <= 1:
        raise ValueError('lmbda must be in range 0-1')
    if not ti_td >= 0:
        raise ValueError('ti_td must be positive')

    new_sequence = deepcopy(sequence)
    for locus in new_sequence.loci:
        models.make_subs_in_locus(locus, q, t=t)
        for i in range(len(locus.codons)):
            if uniform(0, 1) <= lmbda:
                models.make_indel(locus, index=i, ti_td=ti_td, codon_freq=indel_codon_freq)
    return new_sequence

def evolve_tree(sequence,
    taxa=10,
    t=0.01, 
    omega=1.0, 
    kappa=2.0, 
    lmbda=0.001,
    ti_td=0.1,
    codon_freq=None, 
    scale_q=True, 
    model='simple_goldman',
    log=False,
    ):

    """
    Evolve a parent DNA sequence into a set of daughter sequences (taxa) by:
        1. generating a random phylogenetic tree
        2. intantiating a mutational model (e.g. Goldman-Yang-like by default) represented by Q-matrix, with indels
        3. mutate sequence according to tree shape using model 

    Args:
        sequence: a model.Sequence instance
        taxa: number of daughter sequences to evolve
        t: evolution time or branch length
        omega: dN/dS 
        kappa: ratio of transition to transversion rates
        lmbda: probability of indel at codon
        ti_td: ratio of insertions to deletions
        codon_freq: dictionary of codon_frequencies, also know as equilibrium
        frequencies
        scale_q: scales Q so that the average rate of substitution at
        equilibrium equals 1. Branch lengths are thus expected number of nucleotide
        substitutions per codon. See Goldman (1994).
        model: mutational model, 'simple_goldman' will use a Goldman-Yang-like model

    Returns:
        tree instance populated with new sequence strings
    """
    qfunc = model_qfuncs[model]
    q = qfunc(kappa=kappa, omega=omega, codon_freq=codon_freq, scale_q=scale_q, return_dict=False)
   
    tree = trees.random_tree(taxa)
    tree.value = sequence
    for node in trees.get_list_of_tree_nodes(tree)[1:]:
        node.value = evolve_sequence_with_q(node.parent.value, q, t=t, lmbda=lmbda, ti_td=ti_td)
    return tree 
        
def evolve(sequence,
    taxa=10,
    t=0.01, 
    omega=1.0, 
    kappa=2.0, 
    lmbda=0.001,
    ti_td=0.1,
    codon_freq=None, 
    scale_q=True, 
    model='simple_goldman',
    log=False,
    ):
    """
    Wrapper around evolve_tree(). Returns a list of evolved sequences.

    Args:
        sequence: string of DNA nucleotides
        taxa: number of daughter sequences to evolve
        t: evolution time or branch length
        omega: dN/dS 
        kappa: ratio of transition to transversion rates
        lmbda: probability of indel at codon
        ti_td: ratio of insertions to deletions
        codon_freq: dictionary of codon_frequencies, also know as equilibrium
        frequencies
        scale_q: scales Q so that the average rate of substitution at
        equilibrium equals 1. Branch lengths are thus expected number of nucleotide
        substitutions per codon. See Goldman (1994).
        model: mutational model, 'simple_goldman' will use a Goldman-Yang-like model
    """
    sequence = models.Sequence(seq=sequence.lower()) 
    tree = evolve_tree(**locals())
    return [i.value.seq for i in trees.get_list_of_tree_leaves(tree)]


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

def compile_histories(tree):
    """
    Args:
        tree: Tree instance with models.Sequence as value
    """
    if not isinstance(tree.value, models.Sequence):
        raise ValueError('tree.value must be instance of models.Sequence')
    histories = []
    nodes = trees.get_list_of_tree_nodes(tree)
    for t in nodes:
        history = [t.id, t.value.history]
        histories.append(history)
    return histories
 

if __name__=='__main__':
    pass 
