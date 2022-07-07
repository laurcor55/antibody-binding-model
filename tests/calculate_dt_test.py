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
    dock_rotations = [[0, 0, 0]]

    self.ligand_1 = mol.Ligand([50, 50, 50], 9, [0, 0, 0], dock_rotations)
    self.substrates_1 = [mol.Substrate([50, 50, 80], 7, [0, 0.1*np.pi, 0.5*np.pi], dock_rotations)]
    self.reaction_1 = reac.Reaction(self.ligand_1, self.substrates_1, rend, minimum_binding_docks)
    
    self.ligand_2 = mol.Ligand([74, 50, 60], 9, [0, 0, 0], dock_rotations)
    self.substrates_2 = [mol.Substrate([50, 50, 80], 7, [0, 0.1*np.pi, 0.5*np.pi], dock_rotations)]
    self.reaction_2 = reac.Reaction(self.ligand_2, self.substrates_2, rend, minimum_binding_docks)

    self.ligand_3 = mol.Ligand([74, 50, 60], 9, [0, 0, 0], dock_rotations)
    self.substrates_3 = [mol.Substrate([50, 50, 180], 7, [0, 0.1*np.pi, 0.5*np.pi], dock_rotations), mol.Substrate([50, 50, 80], 7, [0, 0.1*np.pi, 0.5*np.pi], dock_rotations)]
    self.reaction_3 = reac.Reaction(self.ligand_3, self.substrates_3, rend, minimum_binding_docks)

  def test_calculate_dt(self):
    close_dock = 3
    substrate_ind = 0
    dt = manual_calculate_dt(self.ligand_1, self.substrates_1, close_dock, substrate_ind)
    self.reaction_1.calculate_dt()
    self.assertAlmostEqual(self.reaction_1.dt,  dt, 10,'not_right')
    
    close_dock = 2
    substrate_ind = 0
    dt = manual_calculate_dt(self.ligand_2, self.substrates_2, close_dock, substrate_ind)
    self.reaction_2.calculate_dt()
    self.assertAlmostEqual(self.reaction_2.dt,  dt, 10,'not_right')

    close_dock = 2
    substrate_ind = 1
    dt = manual_calculate_dt(self.ligand_3, self.substrates_3, close_dock, substrate_ind)
    self.reaction_3.calculate_dt()
    self.assertAlmostEqual(self.reaction_3.dt,  dt, 10,'not_right')

def suite():
    suite = unittest.TestSuite()
    suite.addTest(ReactionTestCase('test_calculate_dt'))
    return suite
    
def manual_calculate_dt(ligand, substrates, close_dock, substrate_ind):
  distance = (np.linalg.norm(ligand.dock_locations[close_dock] - substrates[substrate_ind].dock_locations[close_dock]))*1e-10 
  temperature = 298
  boltzmann = 1.380649e-23
  viscocity = 8.9e-4
  D_0 = boltzmann * temperature /(6*np.pi*viscocity*ligand.radius*1e-10) #m2/s
  D_1 = boltzmann * temperature /(6*np.pi*viscocity*substrates[substrate_ind].radius*1e-10) #m2/s
  D = D_0 + D_1
  return 1/(12*D)*(distance/5)**2

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())