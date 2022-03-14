import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

def get_R(theta, phi, eta):
  R_theta = np.array([[np.cos(theta), -np.sin(theta), 0], \
    [np.sin(theta), np.cos(theta), 0], \
    [0, 0, 1]])

  R_phi = np.array([[np.cos(phi), 0, np.sin(phi)], \
    [0, 1, 0], \
    [-np.sin(phi), 0, np.cos(phi)]])

  R_eta = np.array([[1, 0, 0], \
    [0, np.cos(eta), -np.sin(eta)], \
    [0, np.sin(eta), np.cos(eta)]])    
  
  return np.matmul(np.matmul(R_eta, R_phi), R_theta)


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
    self.ax.set_xlim(-5, 5)
    self.ax.set_ylim(-5, 5)
    self.ax.set_zlim(-5, 5)
    self.update_plot()
  
  def update_plot(self):
    for molecule in self.molecules:
      self.ax.scatter3D(molecule.location[0], molecule.location[1], molecule.location[2])
      self.ax.scatter3D(molecule.dock_offset[0], molecule.dock_offset[1], molecule.dock_offset[2])

    plt.draw()
    plt.pause(1)
  
  def move_molecules(self):
    for molecule in self.molecules:
     # molecule.move()
      molecule.rotate()


class Ligand:
  def __init__(self):
    self.location = np.array([0, 0, 0])
    self.R = get_R(np.random.uniform()*2*np.pi, np.random.uniform()*np.pi, np.random.uniform()*2*np.pi)
    self.dock_offset = np.matmul(np.array([0, 0, 1]), self.R) 

  
  def move(self):
    D = 136e6 # nm2/s
    dt = 10*1e-9 # s
    S = math.sqrt(2*D*dt)
    self.location['x'] += S*np.random.normal()
    self.location['y'] += S*np.random.normal()
    self.location['z'] += S*np.random.normal()
    # Need to multiply S * np.random.normal by R

  def rotate(self):
    D_r = 3.16e7 #/s
    dt = 10*1e-9 # s
    S = math.sqrt(2*D_r*dt)

    delta_eta = S * np.random.normal()
    delta_phi = S * np.random.normal()
    delta_theta = S * np.random.normal()

    R = get_R(delta_theta, delta_phi, delta_eta)
    self.R = np.matmul(self.R, R)
    self.dock_offset = np.matmul(self.R, self.dock_offset)
    print('---')
    print(self.location)
    print(self.dock_offset)
    print(np.linalg.norm(self.dock_offset))

    

molecules = [Ligand()]
Reaction(molecules)