import sys
sys.path.insert(1, '../')
import reaction as reac
import molecules as mol
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import matplotlib.pyplot as plt
import unittest
import matplotlib.pyplot as plt

class ReactionTestCase(unittest.TestCase):
  def setUp(self):
    self.molecules_1 = [mol.Substrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 66, 66], 7, [0, np.pi, -np.pi/2])]

    self.molecules_2 = [mol.FixedSubstrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 66, 66], 7, [0, np.pi, -np.pi/2])]

    self.molecules_3 = [mol.Ligand([50, 66, 66], 7, [0, np.pi, -np.pi/2]), mol.Substrate([50, 50, 50], 9, [0, 0, 0])]

    self.molecules_4 = [mol.Ligand([50, 66, 66], 7, [0, np.pi, -np.pi/2]), mol.Ligand([90, 66, 66], 7, [0, np.pi, -np.pi/2]), mol.Substrate([50, 50, 50], 9, [0, 0, 0]), mol.Substrate([90, 50, 50], 9, [0, 0, 0])]

    self.molecules_5 = [mol.Ligand([50, 66, 66], 7, [0, np.pi, -np.pi/2]), mol.Ligand([90, 66, 66], 7, [0, np.pi, -np.pi/2]), mol.Substrate([50, 50, 50], 9, [0, 0, 0]), mol.Substrate([0, 50, 50], 9, [0, 0, 0])]
    
    self.molecules_6 = [mol.Ligand([50, 66, 2], 7, [0, np.pi, -np.pi/2]), mol.Ligand([90, 66, 66], 7, [0, np.pi, -np.pi/2]), mol.Substrate([50, 50, 50], 9, [0, 0, 0]), mol.Substrate([0, 50, 50], 9, [0, 0, 0])]


  def test_get_locked_pairs(self):
    rend = 1000
    minimum_binding_docks = 3

    self.reaction = reac.Reaction(self.molecules_1, rend, minimum_binding_docks)
    locked_pairs = self.reaction.get_locked_pairs()
    self.assertEqual(locked_pairs[0][0], self.molecules_1[0], 'bad substrate')
    self.assertEqual(locked_pairs[0][1], self.molecules_1[1], 'bad ligand')

    self.reaction = reac.Reaction(self.molecules_2, rend, minimum_binding_docks)
    locked_pairs = self.reaction.get_locked_pairs()
    self.assertEqual(locked_pairs[0][0], self.molecules_2[0], 'bad substrate')
    self.assertEqual(locked_pairs[0][1], self.molecules_2[1], 'bad ligand')\
    
    self.reaction = reac.Reaction(self.molecules_3, rend, minimum_binding_docks)
    locked_pairs = self.reaction.get_locked_pairs()
    self.assertEqual(locked_pairs[0][0], self.molecules_3[1], 'bad substrate')
    self.assertEqual(locked_pairs[0][1], self.molecules_3[0], 'bad ligand')

    self.reaction = reac.Reaction(self.molecules_4, rend, minimum_binding_docks)
    locked_pairs = self.reaction.get_locked_pairs()
    self.assertEqual(locked_pairs[0][0], self.molecules_4[2], 'bad substrate')
    self.assertEqual(locked_pairs[0][1], self.molecules_4[0], 'bad ligand')
    self.assertEqual(locked_pairs[1][0], self.molecules_4[3], 'bad substrate')
    self.assertEqual(locked_pairs[1][1], self.molecules_4[1], 'bad ligand')

    self.reaction = reac.Reaction(self.molecules_5, rend, minimum_binding_docks)
    locked_pairs = self.reaction.get_locked_pairs()
    self.assertEqual(locked_pairs[0][0], self.molecules_5[2], 'bad substrate')
    self.assertEqual(locked_pairs[0][1], self.molecules_5[0], 'bad ligand')

    self.reaction = reac.Reaction(self.molecules_6, rend, minimum_binding_docks)
    locked_pairs = self.reaction.get_locked_pairs()
    self.assertEqual(locked_pairs, [], 'bad substrate')

    
def suite():
  suite = unittest.TestSuite()
  suite.addTest(ReactionTestCase('test_get_locked_pairs'))
  return suite

if __name__ == '__main__':
  runner = unittest.TextTestRunner()
  runner.run(suite())