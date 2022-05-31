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

    self.molecules_1 = [mol.Ligand([50, 50, 50], 9, [0, 0, np.pi]), mol.Substrate([50, 50, 66], 7, [np.pi/2, np.pi, 0])]
    self.molecules_2 = [mol.Ligand([50, 50, 50], 9, [0, 0, 0, 0, 0, 0]), mol.FixedSubstrate([50, 50, 66], 7, [np.pi, 0, np.pi, 0, 0, 0])]

    self.molecules_5 = [mol.Substrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 66, 66], 7, [0, np.pi, 0])]
    self.molecules_6 = [mol.Substrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 66, 66], 7, [-np.pi/2, np.pi, np.pi/2])]
    self.molecules_7 = [mol.FixedSubstrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 66, 66], 7, [-np.pi/2, np.pi, np.pi/2])]
    self.molecules_8 = [mol.FixedSubstrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 66, 1066], 7, [0, np.pi, -np.pi/2])]
    self.molecules_9 = [mol.FixedSubstrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 170, 1066], 7, [0, np.pi, -np.pi/2])]


  def test_check_molecule_status(self):
    rend = 1000
    minimum_binding_docks = 3

    self.reaction = reac.Reaction(self.molecules_1, rend, minimum_binding_docks, 'end_at_binding')
    self.reaction.check_molecule_status()
    self.assertEqual(self.reaction.reaction_status, 'binding', 'reaction bound early')
    self.assertEqual(self.reaction.molecules[0].binding_partner, 2, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[1].binding_partner, 1, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[0].locked_partner, 2, 'wrong locked_parter')
    self.assertEqual(self.reaction.molecules[1].locked_partner, 1, 'wrong locked_parter')
  
    self.reaction = reac.Reaction(self.molecules_2, rend, minimum_binding_docks, 'end_at_binding')
    self.reaction.check_molecule_status()
    self.assertEqual(self.reaction.reaction_status, 'binding', 'reaction status wrong')
    self.assertEqual(self.reaction.molecules[0].binding_partner, 2, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[1].binding_partner, 1, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[0].locked_partner, 2, 'wrong locked_parter')
    self.assertEqual(self.reaction.molecules[1].locked_partner, 1, 'wrong locked_parter')

    self.reaction = reac.Reaction(self.molecules_5, rend, minimum_binding_docks, 'end_at_binding')
    self.reaction.check_molecule_status()
    self.assertEqual(self.reaction.reaction_status, 'free', 'reaction bound early')
    self.assertEqual(self.reaction.molecules[0].binding_partner, 0, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[1].binding_partner, 0, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[0].locked_partner, 0, 'wrong locked_parter')
    self.assertEqual(self.reaction.molecules[1].locked_partner, 0, 'wrong locked_parter')

    self.reaction = reac.Reaction(self.molecules_6, rend, minimum_binding_docks, 'end_at_binding')
    self.reaction.check_molecule_status()
    self.assertEqual(self.reaction.reaction_status, 'free', 'wrong reaction_status')
    self.assertEqual(self.reaction.molecules[0].binding_partner, 0, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[1].binding_partner, 0, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[0].locked_partner, 2, 'wrong locked_parter')
    self.assertEqual(self.reaction.molecules[1].locked_partner, 1, 'wrong locked_parter')

    self.reaction = reac.Reaction(self.molecules_7, rend, minimum_binding_docks, 'end_at_binding')
    self.reaction.check_molecule_status()
    self.assertEqual(self.reaction.reaction_status, 'free', 'wrong reaction_status')
    self.assertEqual(self.reaction.molecules[0].binding_partner, 0, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[1].binding_partner, 0, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[0].locked_partner, 2, 'wrong locked_parter')
    self.assertEqual(self.reaction.molecules[1].locked_partner, 1, 'wrong locked_parter')

    minimum_binding_docks = 2
    self.reaction = reac.Reaction(self.molecules_7, rend, minimum_binding_docks, 'end_at_binding')
    self.reaction.check_molecule_status()
    self.assertEqual(self.reaction.reaction_status, 'binding', 'reaction bound early')
    self.assertEqual(self.reaction.molecules[0].binding_partner, 2, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[1].binding_partner, 1, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[0].locked_partner, 2, 'wrong locked_parter')
    self.assertEqual(self.reaction.molecules[1].locked_partner, 1, 'wrong locked_parter')

    self.reaction = reac.Reaction(self.molecules_8, rend, minimum_binding_docks, 'end_at_binding')
    self.reaction.check_molecule_status()
    self.assertEqual(self.reaction.reaction_status, 'free', 'reaction bound early')
    self.assertEqual(self.reaction.molecules[0].binding_partner, 0, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[1].binding_partner, 0, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[0].locked_partner, 0, 'wrong locked_parter')
    self.assertEqual(self.reaction.molecules[1].locked_partner, 0, 'wrong locked_parter')

    self.reaction = reac.Reaction(self.molecules_9, rend, minimum_binding_docks, 'end_at_binding')
    self.reaction.check_molecule_status()
    self.assertEqual(self.reaction.reaction_status, False, 'reaction bound early')
    self.assertEqual(self.reaction.molecules[0].binding_partner, 0, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[1].binding_partner, 0, 'wrong binding_parter')
    self.assertEqual(self.reaction.molecules[0].locked_partner, 0, 'wrong locked_parter')
    self.assertEqual(self.reaction.molecules[1].locked_partner, 0, 'wrong locked_parter')

    
def suite():
  suite = unittest.TestSuite()
  suite.addTest(ReactionTestCase('test_check_molecule_status'))
  return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())