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
    radius = 18
    ligand_location = [0, 40, 37]
    substrate_location_2 = [0, 50, 0]

    ligand_orientation = [0, np.pi, np.pi/6]
    substrate_orientation = [0, 0, 0]
    dock_rotations = [[0, 0, 0], [0, np.pi, 0]]

    self.ligand_1 = mol.Ligand(ligand_location, radius, ligand_orientation, dock_rotations)
    self.substrate_1 = mol.FixedSubstrate(substrate_location_2, radius, substrate_orientation, dock_rotations)

  def test_find_locked_location(self):
    rend = 1000
    minimum_binding_docks = 3

    self.reaction = reac.Reaction(self.ligand_1, [self.substrate_1], rend, minimum_binding_docks)
    self.reaction.find_locked_location(self.substrate_1, self.ligand_1)
    print(self.reaction.substrates[0].new_location)
    print(self.reaction.ligand.new_location)
    print(self.reaction.ligand.location)

    if self.reaction.ligand.new_location == self.reaction.ligand.location:
      print('hi')


def suite():
  suite = unittest.TestSuite()
  suite.addTest(ReactionTestCase('test_find_locked_location'))
  return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
