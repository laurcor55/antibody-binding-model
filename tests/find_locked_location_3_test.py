import sys
sys.path.insert(1, '../')
import reaction as reac
import molecules as mol
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import matplotlib.pyplot as plt
import unittest
import matplotlib.pyplot as plt
import time

class ReactionTestCase(unittest.TestCase):
  def setUp(self):
    print('hi')

  def test_move_locked_molecules(self):
    rend = 1000
    minimum_binding_docks = 3
    closer_count = 0
    further_count = 0
    seed = 5
    np.random.seed(seed)
    for ii in range(10000):
      self.molecules = [mol.Substrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 66, 66], 7, [0, np.pi, -np.pi/2])]
      molecules = [mol.Substrate([50, 50, 50], radius, [0, 0, 0]), mol.Ligand([50, 43, 87], radius, [0, np.pi, np.pi/8])]
      self.reaction = reac.Reaction(self.molecules, rend, minimum_binding_docks)
      self.reaction.find_locked_location_3(self.reaction.molecules[0], self.reaction.molecules[1])
      old_distances = reac.calculate_distance_docks(self.reaction.molecules[0].dock_locations, self.reaction.molecules[1].dock_locations)
      new_distances = reac.calculate_distance_docks(self.reaction.molecules[0].new_dock_locations, self.reaction.molecules[1].new_dock_locations)
      difference_moved = np.sum(np.array(new_distances)[0:1]) - np.sum(np.array(old_distances)[0:1])
      if difference_moved > 0:
        further_count += 1
      else:
        closer_count += 1
      del self.reaction
      del self.molecules
    print(closer_count)
    print(further_count)

    
def suite():
  suite = unittest.TestSuite()
  suite.addTest(ReactionTestCase('test_move_locked_molecules'))
  return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())