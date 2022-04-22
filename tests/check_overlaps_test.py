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
    self.molecules_1 = [mol.Ligand([50, 50, 50], 9, [0, 0, 0]), mol.Substrate([50, 50, 66], 7, [np.pi, 0, np.pi])]
    self.molecules_2 = [mol.Ligand([50, 50, 50], 9, [0, 0, 0]), mol.FixedSubstrate([50, 50, 66], 7, [np.pi, 0, np.pi])]
    self.molecules_3 = [mol.Substrate([50, 50, 50], 9, [0, 0, 0]), mol.Substrate([50, 50, 66], 7, [np.pi, 0, np.pi])]
    self.molecules_4 = [mol.Ligand([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 50, 66], 7, [np.pi, 0, np.pi])]
    self.molecules_5 = [mol.Substrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 66, 66], 7, [0, np.pi, 0])]
    self.molecules_6 = [mol.Substrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 66, 66], 7, [0, np.pi, -np.pi/2])]
    self.molecules_7 = [mol.FixedSubstrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 66, 66], 7, [0, np.pi, -np.pi/2])]
    self.molecules_8 = [mol.FixedSubstrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 66, 1066], 7, [0, np.pi, -np.pi/2])]
    self.molecules_9 = [mol.FixedSubstrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 170, 1066], 7, [0, np.pi, -np.pi/2])]
    self.molecules_10 = [mol.FixedSubstrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 55, 55], 7, [0, np.pi, -np.pi/2])]
    self.molecules_11 = [mol.Substrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 55, 55], 7, [0, np.pi, -np.pi/2])]
    self.molecules_12 = [mol.Ligand([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 55, 55], 7, [0, np.pi, -np.pi/2])]
    self.molecules_13 = [mol.Substrate([50, 50, 50], 9, [0, 0, 0]), mol.Substrate([50, 55, 55], 7, [0, np.pi, -np.pi/2])]



  def run_check_overlaps(self, molecules):
    rend = 1000
    minimum_binding_docks = 3

    self.reaction = reac.Reaction(molecules, rend, minimum_binding_docks)
    self.reaction.molecules[0].new_location = self.reaction.molecules[0].location
    self.reaction.molecules[1].new_location = self.reaction.molecules[1].location
    return self.reaction.check_overlaps()


  def test_check_reaction_status(self):
    overlaps = self.run_check_overlaps(self.molecules_1)
    self.assertEqual(overlaps, False, 'overlaps wrong')
    overlaps = self.run_check_overlaps(self.molecules_2)
    self.assertEqual(overlaps, False, 'overlaps wrong')
    overlaps = self.run_check_overlaps(self.molecules_3)
    self.assertEqual(overlaps, False, 'overlaps wrong')
    overlaps = self.run_check_overlaps(self.molecules_4)
    self.assertEqual(overlaps, False, 'overlaps wrong')
    overlaps = self.run_check_overlaps(self.molecules_5)
    self.assertEqual(overlaps, False, 'overlaps wrong')
    overlaps = self.run_check_overlaps(self.molecules_6)
    self.assertEqual(overlaps, False, 'overlaps wrong')
    overlaps = self.run_check_overlaps(self.molecules_7)
    self.assertEqual(overlaps, False, 'overlaps wrong')
    overlaps = self.run_check_overlaps(self.molecules_8)
    self.assertEqual(overlaps, False, 'overlaps wrong')
    overlaps = self.run_check_overlaps(self.molecules_9)
    self.assertEqual(overlaps, False, 'overlaps wrong')
    overlaps = self.run_check_overlaps(self.molecules_10)
    self.assertEqual(overlaps, True, 'overlaps wrong')
    overlaps = self.run_check_overlaps(self.molecules_11)
    self.assertEqual(overlaps, True, 'overlaps wrong')
    overlaps = self.run_check_overlaps(self.molecules_12)
    self.assertEqual(overlaps, True, 'overlaps wrong')
    overlaps = self.run_check_overlaps(self.molecules_13)
    self.assertEqual(overlaps, True, 'overlaps wrong')
    
def suite():
  suite = unittest.TestSuite()
  suite.addTest(ReactionTestCase('test_check_reaction_status'))
  return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())