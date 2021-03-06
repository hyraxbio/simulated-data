import trees
import models
import codon_frequencies
from random import uniform, seed
from copy import deepcopy
seed()


def evolve_sequence_with_q(sequence, q, t=1e-2, lmbda=1e-4, ti_td=0.1, indel_codon_freq=None):
    if not 0 <= lmbda <= 1:
        raise ValueError('lmbda must be in range 0-1')
    if not ti_td >= 0:
        raise ValueError('ti_td must be positive')

    new_sequence = deepcopy(sequence)
    p = models.convert_q_to_p(q, t=t)
    p_cumsum, p_codons, p_cumsum_dict = models.get_cumulative_p(p, return_dict=True)
    for locus in new_sequence.loci:
        models.make_subs_in_locus(locus, p_cumsum_dict)
        for i in range(len(locus.codons)):
            if uniform(0, 1) <= lmbda:
                models.make_indel(locus, index=i, ti_td=ti_td, codon_freq=indel_codon_freq)
    return new_sequence


def evolve_tree(sequence,
                taxa=10,
                t=1e-2,
                omega=1.0, 
                kappa=2.0, 
                lmbda=1e-5,
                ti_td=0.1,
                codon_freq='F1x4', 
                scale_q=True, 
                **kwargs
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
        codon_freq: codon frequency model, also know as equilibrium frequencies (default is F1x4)
        scale_q: scales Q so that the average rate of substitution at
        equilibrium equals 1. Branch lengths are thus expected number of nucleotide
        substitutions per codon. See Goldman (1994).

    Returns:
        tree instance populated with new sequence strings
    """

    codon_freq = getattr(codon_frequencies, codon_freq)(sequence)

    q = models.goldman_Q(kappa=kappa, omega=omega, codon_freq=codon_freq, scale_q=scale_q, return_dict=False)
   
    tree = trees.random_tree(taxa)
    tree.value = sequence
    for node in trees.get_list_of_tree_nodes(tree)[1:]:
        node.value = evolve_sequence_with_q(node.parent.value, q, t=t, lmbda=lmbda, ti_td=ti_td)
    return tree 
       
 
def evolve(sequence,
           taxa=10,
           t=1e-2,
           omega=1.0, 
           kappa=2.0, 
           lmbda=1e-4,
           ti_td=0.1,
           codon_freq='F1x4', 
           scale_q=True, 
           strip_deletions=False,
           log=False,
           verbose=False,
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
        codon_freq: codon frequency model, also know as equilibrium frequencies (default is F1x4)
        scale_q: scales Q so that the average rate of substitution at
        equilibrium equals 1. Branch lengths are thus expected number of nucleotide
        substitutions per codon. See Goldman (1994).
        log: if True, returns list of evolved sequences AND list of mutations
        verbose: if True, prints parameters
        strip_deletions: False,
    """
    if not isinstance(sequence, str):
        raise TypeError('sequence must be a string')
   
    if verbose: 
        print '\n\n---------------------------------------'
        print 'evolving new sequences with parameters:'
        print '---------------------------------------'
        for i, j in zip(['taxa',
                         't',
                         'omega',
                         'kappa',
                         'lmbda',
                         'ti_td',
                         'scale_q',
                         'codon_freq',
                         'strip_deletions',
                         'log',
                         ],
                        [taxa, t, omega, kappa, lmbda, ti_td, scale_q, codon_freq, strip_deletions, log]):
            print '{0:<10} {1:<15}'.format(i, j)
        print '---------------------------------------\n\n'

    sequence = models.Sequence(seq=sequence.lower()) 
    tree = evolve_tree(**locals())
    leaves = trees.get_list_of_tree_leaves(tree)
    sequences = [i.value.seq for i in leaves]
    if strip_deletions:
        sequences = [i.replace('-', '') for i in sequences]
    if log:
        mutations = [i.value.mutations for i in leaves]
        return sequences, mutations
    return sequences


def print_mutations(old_sequence, new_sequence, colour=True):
    """
    Simple string comparison.
    """

    sim = ''
    for i, j in zip(old_sequence, new_sequence):
        if i == j:
            sim += ' '
        else:
            sim += '|'

    old_sequence_aa = [i for i in models.convert_dna_to_aa(old_sequence)] 
    new_sequence_aa = [i for i in models.convert_dna_to_aa(new_sequence)]
  
    if colour:
        red = ['\x1b[31m', '\x1b[0m']
        for i in range(len(old_sequence_aa)):
            if old_sequence_aa[i] != new_sequence_aa[i]:
                new_sequence_aa[i] = red[0] + new_sequence_aa[i] + red[1]
                
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
    histories = {}
    nodes = trees.get_list_of_tree_nodes(tree)
    for t in nodes:
        histories[t.id] = t.value.history
    return histories


if __name__ == '__main__':
    pass
