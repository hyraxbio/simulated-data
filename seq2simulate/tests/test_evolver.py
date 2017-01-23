import unittest
import evolver

class TreeTester(unittest.TestCase):
  
    def test_inits(self):
        tree0 = evolver.Tree()
        self.assertIsInstance(tree0.id, int)
        self.assertIsNone(tree0.value)
        self.assertIsNone(tree0.left)
        self.assertIsNone(tree0.right)
        self.assertEqual(tree0.__repr__(), '<Tree {}: val:None len:0.1 left:None right:None parent:None>'.format(tree0.id))

    def test_tree_ids_increment_automatically(self):
        tree0 = evolver.Tree()
        tree1 = evolver.Tree()
        self.assertEqual(tree1.id, tree0.id + 1)

    def test_get_two_random_orphans(self):
        nodes = [evolver.Tree() for i in range(4)]
        nodes[0].left = nodes[1]
        nodes[0].right = nodes[2]
        self.assertEqual(nodes[1].parent, nodes[0])
        self.assertEqual(nodes[2].parent, nodes[0])
        orphans = evolver.get_two_random_orphans(nodes)
        self.assertEqual(len(orphans), 2)
        for i in orphans:
            self.assertIsNone(i.parent) 

    def test_get_orphans(self):
        nodes = [evolver.Tree() for i in range(4)]
        self.assertEqual(len(evolver.get_orphans(nodes)), 4)
        nodes[0].parent = nodes[1]
        self.assertEqual(len(evolver.get_orphans(nodes)), 3)

    def test_preorder_exec(self):
        t = evolver.random_tree(10)  
        tlist = evolver.get_list_of_tree_nodes(t)
        self.assertTrue(all(i.get_value() is None for i in tlist))
        evolver.preorder_exec(t, function='set_value', arguments=['test'])
        self.assertTrue(all(i.get_value() == 'test' for i in tlist))
 
    def test_preorder_exec_validation(self):
        t = evolver.random_tree(10)  
        with self.assertRaises(ValueError):
            evolver.preorder_exec('blah')
        with self.assertRaises(ValueError):
            evolver.preorder_exec(t, function=5)
        with self.assertRaises(ValueError):
            evolver.preorder_exec(t, arguments=5)
        
    def test_get_list_of_tree_nodes(self):
        t = evolver.random_tree(10)  
        tl = evolver.get_list_of_tree_nodes(t)
        self.assertEqual(len(tl), 19)

    def test_get_list_of_tree_nodes_validation(self):
        t = evolver.random_tree(10)  
        with self.assertRaises(ValueError):
            evolver.get_list_of_tree_nodes('blah')

if __name__=='__main__':
    unittest.main()
