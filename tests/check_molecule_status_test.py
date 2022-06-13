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

    ligand_location = [0, 0, 37]
    substrate_location = [0, 0, 0]
    substrate_location_2 = [0, 50, 0]

    ligand_orientation = [0, np.pi, 0]
    substrate_orientation = [0, 0, 0]
    n_docks = 4
    self.ligand_1 = mol.Ligand(ligand_location, radius, ligand_orientation, n_docks)
    self.substrates_1 = [mol.FixedSubstrate(substrate_location, radius, substrate_orientation, n_docks), mol.FixedSubstrate(substrate_location_2, radius, substrate_orientation, n_docks)]
    
  def test_check_molecule_status(self):
    rend = 1000
    minimum_binding_docks = 3
    reaction = reac.Reaction(self.ligand_1, self.substrates_1, rend, minimum_binding_docks)
    reaction.check_molecule_status()
    self.assertEqual(reaction.reaction_status, 'binding', 'not binding')
    self.assertEqual(reaction.substrates[0].locked, True, 'not locked')
    self.assertEqual(reaction.ligand.locked, True, 'not locked')


def suite():
  suite = unittest.TestSuite()
  suite.addTest(ReactionTestCase('test_check_molecule_status'))
  return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())