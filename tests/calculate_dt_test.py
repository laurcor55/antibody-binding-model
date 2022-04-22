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
    rstart = 30
    rend = 1000
    minimum_binding_docks = 3

    self.molecule_a = mol.Ligand(np.array([50, 50, 50], float), 9, [0, 0, 0])
    self.molecule_b = mol.Substrate(np.array([50, 50, 80], float), 7, [0, 0.1*np.pi, 0.5*np.pi])
    self.molecule_c = mol.Ligand(np.array([75, 50, 60], float), 9, [0, 0, 0])

    self.molecules_1 = [self.molecule_a, self.molecule_b]
    self.reaction_1 = reac.Reaction(self.molecules_1, rend, minimum_binding_docks)
    
    self.molecules_2 = [self.molecule_b, self.molecule_c]
    self.reaction_2 = reac.Reaction(self.molecules_2, rend, minimum_binding_docks)

    self.molecules_3 = [self.molecule_a, self.molecule_b, self.molecule_c]
    self.reaction_3 = reac.Reaction(self.molecules_3, rend, minimum_binding_docks)
  
  def test_calculate_dt(self):
    close_dock = 1
    close_molecule_1 = 0
    close_molecule_2 = 1
    dt = manual_calculate_dt(self.molecules_1, close_dock, close_molecule_1, close_molecule_2)
    self.reaction_1.calculate_dt()
    self.assertAlmostEqual(self.reaction_1.dt,  dt, 'incorrect calculate_dt')

    close_dock = 0
    close_molecule_1 = 0
    close_molecule_2 = 1
    dt = manual_calculate_dt(self.molecules_2, close_dock, close_molecule_1, close_molecule_2)
    self.reaction_2.calculate_dt()
    self.assertAlmostEqual(self.reaction_2.dt,  dt, 'incorrect calculate_dt')

    close_dock = 0
    close_molecule_1 = 1
    close_molecule_2 = 2
    dt = manual_calculate_dt(self.molecules_3, close_dock, close_molecule_1, close_molecule_2)
    self.reaction_3.calculate_dt()
    self.assertAlmostEqual(self.reaction_3.dt,  dt, 'incorrect calculate_dt')

def suite():
    suite = unittest.TestSuite()
    suite.addTest(ReactionTestCase('test_calculate_dt'))
    return suite
    
def manual_calculate_dt(molecules, close_dock, close_molecule_1, close_molecule_2):
  distance = np.linalg.norm(molecules[close_molecule_1].dock_locations[close_dock] - molecules[close_molecule_2].dock_locations[close_dock])*1e-10
  temperature = 298
  boltzmann = 1.380649e-23
  viscocity = 8.9e-4
  D_0 = boltzmann * temperature /(6*np.pi*viscocity*molecules[0].radius*1e-10) #m2/s
  D_1 = boltzmann * temperature /(6*np.pi*viscocity*molecules[1].radius*1e-10) #m2/s
  D = D_0 + D_1
  return 1/(12*D)*(distance/10)**2

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())