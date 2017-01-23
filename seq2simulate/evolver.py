
import random

class Tree(object):
    id = 0

    def __init__(self, value=None, length=0.1, left=None, right=None, parent=None):
        self.id = Tree.id
        Tree.id += 1
        self.value = value
        self.length = length
        self.left  = left
        self.right = right
        self.parent = parent

    @property
    def left(self):
        return self._left

    @left.setter
    def left(self, left):
        self._left = left
        if left is not None:
            left.parent = self

    @property
    def right(self):
        return self._right

    @right.setter
    def right(self, right):
        self._right = right
        if right is not None:
            right.parent = self

    def set_value(self, value):
        self.value = value

    def get_value(self):
        return self.value
    
    def __repr__(self):
        format_str = '<Tree {}: val:{} len:{} left:{} right:{} parent:{}>'
        return format_str.format(*(str(i) for i in [self.id, self.value, self.length, 
            self._get_id_or_none(self.left),
            self._get_id_or_none(self.right),
            self._get_id_or_none(self.parent),
            ]))

    def _get_id_or_none(self, obj):
        try:
            val = getattr(obj, 'id')
        except:
            val = None
        return val
 
def random_tree(num_taxa, mean_branch_length=0.1):
    """
    Returns a randomly generated phylogenetic tree.

    Args:
        num_taxa: The number of taxa.
        mean_branch_length: The average substitution rate of the branches.
            
    """
    nodes = [Tree() for i in range(num_taxa)]
    
    while len(get_orphans(nodes)) > 1: 
        orphan_nodes = get_two_random_orphans(nodes)
        new_node = Tree(left=orphan_nodes[0], right=orphan_nodes[1])
        nodes.append(new_node)

    return get_orphans(nodes)[0]

def get_orphans(tree_nodes):
    """
    Returns orphaned Tree nodes from list of Trees.
    """
    return [i for i in tree_nodes if i.parent is None]

def get_two_random_orphans(tree_nodes):
    """
    Returns two random orphaned Tree nodes from list of Trees.
    """
    orphan_nodes = []
    while len(orphan_nodes) < 2:
        orphan_node = random.choice(tree_nodes)
        if orphan_node in orphan_nodes:
            continue
        if orphan_node.parent is not None:
            continue
        orphan_nodes.append(orphan_node)
    return orphan_nodes

def preorder_exec(tree, function='__str__', arguments=[]):
    """
    Perform a preorder traversal of 'tree', and execute 'function' passing 'arguments'.
    """
    if tree:
        print(tree.__repr__())
        getattr(tree, function)(*arguments)
        preorder_exec(tree.left, function=function, arguments=arguments)
        preorder_exec(tree.right, function=function, arguments=arguments)

def get_list_of_tree_nodes(tree, nodes=[]):
    """
    Perform a preorder traversal of 'tree', and return a list of all nodes.
    """
    if tree:
        nodes.append(tree)
        get_list_of_tree_nodes(tree.left, nodes=nodes)
        get_list_of_tree_nodes(tree.right, nodes=nodes)

    return nodes

if __name__=='__main__':
    pass
    
