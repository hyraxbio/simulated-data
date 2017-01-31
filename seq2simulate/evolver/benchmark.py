import evolve, models, trees
import time

with open('../../data/split/Seq2_Sus', 'r') as f:
    old_sequence = f.read()[10:].replace('\n', '').lower()

old_sequence = models.Sequence(old_sequence)

for i in [3, 6, 9]:
    start_time = time.time()
    #new_sequences, mutations = evolve.evolve(old_sequence, taxa=10, log=True)
    tree = evolve.evolve_tree(old_sequence, taxa=i)
    nodes = len(trees.get_list_of_tree_nodes(tree))
    print("taxa:{} (nodes={}), sec:{}".format(i, nodes, time.time() - start_time))
