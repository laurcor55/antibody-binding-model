import math
import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

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

def update_lines(num, dataLines, lines):
  for line, data in zip(lines, dataLines):
    # NOTE: there is no .set_data() for 3 dim data...
    x = np.array([0, data[0:2, :num]])
    print(x.shape)
    print(data[0:2, :num].shape)
    line.set_data(data[0:2, :num])
    line.set_3d_properties(data[2, :num])
  return lines



class Reaction:
  def __init__(self, molecules, t_steps):
    self.molecules = molecules
    self.t_steps = t_steps
    self.time = 0
    self.progress_reaction()
    self.make_animation_2()
  
  def make_animation(self):
    self.fig = plt.figure()
    self.ax = p3.Axes3D(self.fig)
    data = [np.rot90(np.asarray(molecule.dock_offsets)) for molecule in self.molecules]
    #lines = [self.ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]
    lines = [self.ax.plot([], [], [])[0] for _ in data]
    self.ax.set_xlim3d(-5, 5)
    self.ax.set_ylim3d(-5, 5)
    self.ax.set_zlim3d(-5, 5)
    line_ani = animation.FuncAnimation(self.fig, update_lines, self.t_steps, fargs=(data, lines), interval=50, blit=True)
    plt.show()
    
  def progress_reaction(self):
    for ii in range(self.t_steps):
      self.move_molecules()

  def make_animation_2(self):
    self.fig = plt.figure()
    self.ax = p3.Axes3D(self.fig)
    limits = 50
    for t_step in range(self.t_steps):
      self.ax.clear()
      self.ax.set_xlim3d(-limits, limits)
      self.ax.set_ylim3d(-limits, limits)
      self.ax.set_zlim3d(-limits, limits)
      for molecule in self.molecules:
        x = np.array([molecule.locations[t_step][0], molecule.locations[t_step][0] + molecule.dock_offsets[t_step][0]])
        y = np.array([molecule.locations[t_step][1], molecule.locations[t_step][1] + molecule.dock_offsets[t_step][1]])
        z = np.array([molecule.locations[t_step][2], molecule.locations[t_step][2] + molecule.dock_offsets[t_step][2]])
        self.ax.plot3D(x, y, z)
      plt.pause(0.1)
    plt.show()
  
  def move_molecules(self):
    for molecule in self.molecules:
      molecule.move()
      molecule.rotate()

class Ligand:
  def __init__(self):
    self.location = np.zeros(3)
    self.locations = [self.location.copy()]
    self.R = get_R(np.random.uniform()*2*np.pi, np.random.uniform()*np.pi, np.random.uniform()*2*np.pi)
    self.dock_offset = np.matmul(np.array([0, 0, 1]), self.R) 
    self.dock_offsets = [self.dock_offset]
  
  def move(self):
    D = 136e6 # nm2/s
    dt = 10*1e-9 # s
    S = math.sqrt(2*D*dt)
    self.location[0] += S*np.random.normal()
    self.location[1] += S*np.random.normal()
    self.location[2] += S*np.random.normal()
    self.locations.append(self.location.copy())
    print(self.locations)

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
    self.dock_offsets.append(self.dock_offset)


    
t_steps = 50
molecules = [Ligand(), Ligand()]
Reaction(molecules, t_steps)