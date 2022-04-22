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
    print('need to add')

  def test_move_locked_molecules(self):
    print('need to add')
    
def suite():
  suite = unittest.TestSuite()
  suite.addTest(ReactionTestCase('test_move_locked_molecules'))
  return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())