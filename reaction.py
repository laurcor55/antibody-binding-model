import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

class Reaction:
  def __init__(self, molecules):
    self.molecules = molecules
    self.init_plot()
    self.time = 0
    self.progress_reaction()
    
  def progress_reaction(self):
    for ii in range(10):
      self.move_molecules()
      self.update_plot()
  
  def init_plot(self):
    self.fig = plt.figure()
    self.ax = plt.axes(projection='3d')
    self.ax.set_xlim(0, 20)
    self.ax.set_ylim(0, 20)
    self.ax.set_zlim(0, 20)
    self.update_plot()
  
  def update_plot(self):
    for molecule in self.molecules:
      self.ax.scatter3D(molecule.location['x'], molecule.location['y'], molecule.location['z'])
    plt.draw()
    plt.pause(1)
  
  def move_molecules(self):
    for molecule in self.molecules:
      molecule.move()


class Ligand:
  def __init__(self):
    self.location = {'x':10, 'y':10, 'z':10}
    self.orientation = 5# {eta:0, phi:0, theta:0}
  
  def move(self):
    D = 136e6 # nm2/s
    D_r = 3.16e7 #/s
    dt = 10*1e-9 # s
    S = math.sqrt(2*D*dt)
    self.location['x'] += S*np.random.normal()
    self.location['y'] += S*np.random.normal()
    self.location['z'] += S*np.random.normal()
    print(self.location)
    #theta = math.sqrt(2*D_r*dt)



molecules = [Ligand(), Ligand()]
Reaction(molecules)