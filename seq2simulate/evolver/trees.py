import re
import random
from types import NoneType

class Tree(object):
    id = 0

    def __init__(self, 
        value=None, 
        length=0.1, 
        left=None, 
        right=None, 
        parent=None):
        

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
        # this is convenient but it concerns me...
        if left is not None:
            left.parent = self

    @property
    def right(self):
        return self._right

    @right.setter
    def right(self, right):
        self._right = right
        # this is convenient but it concerns me...
        if right is not None:
            right.parent = self

    @property
    def length(self):
        return self._length

    @length.setter
    def length(self, length):
        if length < 0.0 or length > 1.0:
            raise ValueError('length must be in range 0-1')
        self._length = length

    def set_value(self, value):
        self.value = value

    def get_value(self):
        return self.value
    
    def __str__(self):
        format_str = '<Tree {}: len:{} left:{} right:{} parent:{}>'
        return format_str.format(*(str(i) for i in [self.id, self.length, 
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

    def to_list(self):
        _id = self.id
        left = self.left
        if left is not None:
            left = left.to_list()
        right = self.right
        if right is not None:
            right = right.to_list()
        if left is None and right is None:
            return _id
        else:
            return [_id, [left, right]]

def random_tree(num_taxa, mean_branch_length=0.1):
    """
    Returns a randomly generated phylogenetic tree.

    Args:
        num_taxa: The number of taxa.
        mean_branch_length: The average substitution rate of the branches.
            
    """
    if mean_branch_length < 0.0 or mean_branch_length > 1.0:
        raise ValueError('mean_branch_length must be in range 0-1')

    nodes = [Tree() for i in range(int(num_taxa))]
    
    while len(get_orphans(nodes)) > 1: 
        orphan_nodes = get_two_random_orphans(nodes)
        new_node = Tree(left=orphan_nodes[0], 
                        right=orphan_nodes[1], 
                        length=random.uniform(0, 2*mean_branch_length)
                    )
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
    orphan_nodes = [i for i in tree_nodes if i.parent is None]
    if len(orphan_nodes) < 2:
        return []
    orphans = []
    while len(orphans) < 2:
        orphan_node = random.choice(orphan_nodes)
        if orphan_node in orphans:
            continue
        orphans.append(orphan_node)
    return orphans

def preorder_exec(tree, function='__str__', arguments=[]):
    """
    Perform a preorder traversal of 'tree', and execute 'function' passing 'arguments'.
    """

    if tree:
        getattr(tree, function)(*arguments)
        preorder_exec(tree.left, function=function, arguments=arguments)
        preorder_exec(tree.right, function=function, arguments=arguments)

def get_list_of_tree_nodes(tree):
    """
    Perform a preorder traversal of 'tree', and return a list of all nodes.
    """

    def get_nodes(tree, nodes=[]):

        if tree:
            nodes.append(tree)
            get_nodes(tree.left, nodes=nodes)
            get_nodes(tree.right, nodes=nodes)

        return nodes

    nodes = []
    return get_nodes(tree, nodes=nodes)

def get_dict_of_tree_values(tree):
    """
    Return a dictionary of node values indexed by id.
    """
    nodes = get_list_of_tree_nodes(tree)
    return {i.id:i.value for i in nodes}

def get_list_of_tree_leaves(tree):
    """
    Perform a preorder traversal of 'tree', and return a list of all terminal nodes.
    """

    def get_nodes(tree, nodes=[]):

        if tree:
            if tree.left is None and tree.right is None:
                nodes.append(tree)
            get_nodes(tree.left, nodes=nodes)
            get_nodes(tree.right, nodes=nodes)

        return nodes

    nodes = []
    return get_nodes(tree, nodes=nodes)


if __name__=='__main__':
    pass
