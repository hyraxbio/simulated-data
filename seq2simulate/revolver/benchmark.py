import evolve
import models
import trees
import codon_frequencies
import time


def benchmark_evolve_tree():
    with open('../../data/split/Seq2_Sus', 'r') as f:
        old_sequence = f.read()[10:].replace('\n', '').lower()
    
    old_sequence = models.Sequence(old_sequence)
    
    for i in [3, 6, 9]:
        start_time = time.time()
        # new_sequences, mutations = evolve.evolve(old_sequence, taxa=10, log=True)
        tree = evolve.evolve_tree(old_sequence, taxa=i)
        nodes = len(trees.get_list_of_tree_nodes(tree))
        print("taxa:{} (nodes={}), sec:{}".format(i, nodes, time.time() - start_time))


def benchmark_evolve():
    with open('../../data/split/Seq2_Sus', 'r') as f:
        old_sequence = f.read()[10:].replace('\n', '').lower()
    taxa = 10 
    start_time = time.time()
    new_sequences, mutations = evolve.evolve(old_sequence, taxa=taxa, t=0.01, log=True)
    print("taxa:{}, sec:{}".format(taxa, time.time() - start_time))
    return old_sequence, new_sequences, mutations


def benchmark_evolve_big_set():
    with open('../../data/split/Seq2_Sus', 'r') as f:
        old_sequence = f.read()[10:].replace('\n', '').lower()
    taxa = 100
    start_time = time.time()
    new_sequences, mutations = evolve.evolve(old_sequence, taxa=taxa, log=True, omega=0.1)
    print("taxa:{}, sec:{}".format(taxa, time.time() - start_time))
    return old_sequence, new_sequences, mutations


def plot_p():
    with open('../../data/split/Seq2_Sus', 'r') as f:
        old_sequence = f.read()[10:].replace('\n', '').lower()
    seq = models.Sequence(old_sequence)
    codon_freq = codon_frequencies.F1x4(seq)
    q = models.goldman_Q(codon_freq=codon_freq)
    models.plot_p_over_time(q, t=0.1, codon='aaa', logscale=False)
    q = models.goldman_Q()
    models.plot_p_over_time(q, t=0.1, codon='aaa', logscale=False)


if __name__ == '__main__':
    # benchmark_evolve_tree()
    ols, ns, mt = benchmark_evolve()
    ols, ns, mt = benchmark_evolve_big_set()
    evolve.print_mutations(ols, ns[0])
    # plot_p()

