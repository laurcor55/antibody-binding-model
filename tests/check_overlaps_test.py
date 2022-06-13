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
    ligand_location = [0, 0, 37]
    substrate_location = [0, 0, 0]
    substrate_location_2 = [0, 50, 0]

    ligand_orientation = [0, np.pi, 0]
    substrate_orientation = [0, 0, 0]
    n_docks = 4
    self.ligand_1 = mol.Ligand(ligand_location, radius, ligand_orientation, n_docks)
    self.substrates_1 = [mol.FixedSubstrate(substrate_location, radius, substrate_orientation, n_docks), mol.FixedSubstrate(substrate_location_2, radius, substrate_orientation, n_docks)]

    ligand_location = [0, 0, 30]
    substrate_location = [0, 0, 0]
    substrate_location_2 = [0, 50, 0]

    self.ligand_2 = mol.Ligand(ligand_location, radius, ligand_orientation, n_docks)
    self.substrates_2 = [mol.FixedSubstrate(substrate_location, radius, substrate_orientation, n_docks), mol.FixedSubstrate(substrate_location_2, radius, substrate_orientation, n_docks)]

  def run_check_overlaps(self, ligand, substrates):
    rend = 1000
    minimum_binding_docks = 3

    self.reaction = reac.Reaction(ligand, substrates, rend, minimum_binding_docks)
    self.reaction.ligand.new_location = self.reaction.ligand.location
    self.reaction.substrates[0].new_location = self.reaction.substrates[0].location
    self.reaction.substrates[1].new_location = self.reaction.substrates[1].location
    return self.reaction.check_overlaps()


  def test_check_reaction_status(self):
    overlaps = self.run_check_overlaps(self.ligand_1, self.substrates_1)
    self.assertEqual(overlaps, True, 'overlaps wrong')
  
    overlaps = self.run_check_overlaps(self.ligand_2, self.substrates_2)
    self.assertEqual(overlaps, False, 'overlaps wrong')
  
def suite():
  suite = unittest.TestSuite()
  suite.addTest(ReactionTestCase('test_check_reaction_status'))
  return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())