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

  def test_find_locked_location(self):
    rend = 1000
    minimum_binding_docks = 3

    self.reaction = reac.Reaction(self.molecules_1, rend, minimum_binding_docks)
    substrate = self.molecules_1[0]
    ligand = self.molecules_1[1]
    self.reaction.find_locked_location(substrate, ligand)
    self.reaction.check_reaction_status()
    self.assertEqual(self.reaction.end_reaction, True, 'reaction ended early')
    self.assertEqual(self.reaction.binding, True, 'reaction bound early')
    self.assertEqual(self.reaction.molecules[0].binding_partner, 2, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[1].binding_partner, 1, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[0].locked_partner, 2, 'wrong locked_parter')
    self.assertEqual(self.reaction.molecules[1].locked_partner, 1, 'wrong locked_parter')

    self.reaction = reac.Reaction(self.molecules_2, rend, minimum_binding_docks)
    substrate = self.molecules_2[0]
    ligand = self.molecules_2[1]
    self.reaction.find_locked_location(substrate, ligand)
    self.reaction.check_reaction_status()
    self.assertEqual(self.reaction.end_reaction, True, 'reaction ended early')
    self.assertEqual(self.reaction.binding, True, 'reaction bound early')
    self.assertEqual(self.reaction.molecules[0].binding_partner, 2, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[1].binding_partner, 1, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[0].locked_partner, 2, 'wrong locked_parter')
    self.assertEqual(self.reaction.molecules[1].locked_partner, 1, 'wrong locked_parter')
    
def suite():
  suite = unittest.TestSuite()
  suite.addTest(ReactionTestCase('test_find_locked_location'))
  return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
