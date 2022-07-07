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
    dock_rotations = [[0, 0, 0], [0, np.pi, 0]]
    self.ligand_1 = mol.Ligand(ligand_location, radius, ligand_orientation, dock_rotations)
    self.substrates_1 = [mol.FixedSubstrate(substrate_location, radius, substrate_orientation, dock_rotations), mol.FixedSubstrate(substrate_location_2, radius, substrate_orientation, dock_rotations)]

    radius = 18
    ligand_location = [0, 0, 37]
    substrate_location = [0, 0, 0]
    substrate_location_2 = [0, 50, 0]

    ligand_orientation = [0, np.pi-np.pi/4, 0]
    substrate_orientation = [0, 0, 0]
    dock_rotations = [[0, 0, 0], [0, 0, 0]]
    self.ligand_2 = mol.Ligand(ligand_location, radius, ligand_orientation, dock_rotations)
    self.substrates_2 = [mol.FixedSubstrate(substrate_location, radius, substrate_orientation, dock_rotations), mol.FixedSubstrate(substrate_location_2, radius, substrate_orientation, dock_rotations)]

    radius = 18
    ligand_location = [0, 50, 37]
    substrate_location = [0, 0, 0]
    substrate_location_2 = [0, 50, 0]

    ligand_orientation = [0, np.pi, 0]
    substrate_orientation = [0, 0, 0]
    dock_rotations = [[0, 0, 0]]

    self.ligand_3 = mol.Ligand(ligand_location, radius, ligand_orientation, dock_rotations)
    self.substrates_3 = [mol.FixedSubstrate(substrate_location, radius, substrate_orientation, dock_rotations), mol.FixedSubstrate(substrate_location_2, radius, substrate_orientation, dock_rotations)]

    radius = 18
    ligand_location = [0, 50, 37]
    substrate_location = [0, 0, 0]
    substrate_location_2 = [0, 50, 0]

    ligand_orientation = [0, np.pi, 0]
    substrate_orientation = [0, 0, 0]
    dock_rotations = [[0, 0, 0], [0, np.pi/4, 0]]

    self.ligand_4 = mol.Ligand(ligand_location, radius, ligand_orientation, dock_rotations)
    self.substrates_4 = [mol.FixedSubstrate(substrate_location, radius, substrate_orientation, dock_rotations), mol.FixedSubstrate(substrate_location_2, radius, substrate_orientation, dock_rotations)]

    radius = 18
    ligand_location = [0, 40, 37]
    substrate_location = [0, 0, 0]
    substrate_location_2 = [0, 50, 0]

    ligand_orientation = [0, np.pi, np.pi/6]
    substrate_orientation = [0, 0, 0]
    dock_rotations = [[0, 0, 0], [0, np.pi, 0]]

    self.ligand_5 = mol.Ligand(ligand_location, radius, ligand_orientation, dock_rotations)
    self.substrates_5 = [mol.FixedSubstrate(substrate_location, radius, substrate_orientation, dock_rotations), mol.FixedSubstrate(substrate_location_2, radius, substrate_orientation, dock_rotations)]
    
  def test_check_molecule_status(self):
    rend = 1000
    minimum_binding_docks = 3
    reaction = reac.Reaction(self.ligand_1, self.substrates_1, rend, minimum_binding_docks)
    reaction.check_molecule_status()
  
    self.assertEqual(reaction.reaction_status, 'binding', 'not binding')
    self.assertEqual(reaction.substrates[0].locked, True, 'not locked')
    self.assertEqual(reaction.substrates[1].locked, False, 'not locked')

    reaction = reac.Reaction(self.ligand_2, self.substrates_2, rend, minimum_binding_docks)
    reaction.check_molecule_status()
    self.assertEqual(reaction.reaction_status, 'free', 'not binding')
    self.assertEqual(reaction.substrates[0].locked, False, 'not locked')
    self.assertEqual(reaction.substrates[1].locked, False, 'not locked')
    

    reaction = reac.Reaction(self.ligand_3, self.substrates_3, rend, minimum_binding_docks)
    reaction.check_molecule_status()
    self.assertEqual(reaction.reaction_status, 'binding', 'not binding')
    self.assertEqual(reaction.substrates[0].locked, False, 'not locked')
    self.assertEqual(reaction.substrates[1].locked, True, 'not locked')

    reaction = reac.Reaction(self.ligand_4, self.substrates_4, rend, minimum_binding_docks)
    reaction.check_molecule_status()
    self.assertEqual(reaction.reaction_status, 'binding', 'not binding')
    self.assertEqual(reaction.substrates[0].locked, False, 'not locked')
    self.assertEqual(reaction.substrates[1].locked, True, 'not locked')

    reaction = reac.Reaction(self.ligand_5, self.substrates_5, rend, minimum_binding_docks)
    reaction.check_molecule_status()
    self.assertEqual(reaction.reaction_status, 'free', 'not binding')
    self.assertEqual(reaction.substrates[0].locked, False, 'not locked')
    self.assertEqual(reaction.substrates[1].locked, True, 'not locked')


def suite():
  suite = unittest.TestSuite()
  suite.addTest(ReactionTestCase('test_check_molecule_status'))
  return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())