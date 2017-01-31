import evolve, models, trees
import time

def benchmark_evolve_tree():
    with open('../../data/split/Seq2_Sus', 'r') as f:
        old_sequence = f.read()[10:].replace('\n', '').lower()
    
    old_sequence = models.Sequence(old_sequence)
    
    for i in [3, 6, 9]:
        start_time = time.time()
        #new_sequences, mutations = evolve.evolve(old_sequence, taxa=10, log=True)
        tree = evolve.evolve_tree(old_sequence, taxa=i)
        nodes = len(trees.get_list_of_tree_nodes(tree))
        print("taxa:{} (nodes={}), sec:{}".format(i, nodes, time.time() - start_time))

def benchmark_evolve():
    with open('../../data/split/Seq2_Sus', 'r') as f:
        old_sequence = f.read()[10:].replace('\n', '').lower()
    taxa = 10 
    start_time = time.time()
    new_sequences, mutations = evolve.evolve(old_sequence, taxa=taxa, log=True)
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

if __name__=='__main__':
    benchmark_evolve_tree()
    ols, ns, mt = benchmark_evolve()
    ols, ns, mt = benchmark_evolve_big_set()
