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
    self.molecules = [mol.Substrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 66, 66], 7, [0, np.pi, -np.pi/2])]

  def test_move_locked_molecules(self):
    rend = 1000
    minimum_binding_docks = 3

    self.reaction = reac.Reaction(self.molecules, rend, minimum_binding_docks)
    self.reaction.move_locked_molecules()
    
def suite():
  suite = unittest.TestSuite()
  suite.addTest(ReactionTestCase('test_move_locked_molecules'))
  return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())