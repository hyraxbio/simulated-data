import trees, models

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

def print_mutations(old_sequence, new_sequence):
    sim = ''
    for i,j in zip(old_sequence, new_sequence):
        if i==j:
            sim += ' '
        else:
            sim += '|'
    print(old_sequence)
    print(sim)
    print(new_sequence)

if __name__=='__main__':
    #t10 = trees.random_tree(10)
    old_sequence = 'atgcaacggcgattatacgtatcgtgcatcgatcatcgcatgcaacggcgattatacgtatcgtgcatcgatcatcgc'
    new_sequence = evolve_sequence(old_sequence) 
    print_mutations(old_sequence, new_sequence)
