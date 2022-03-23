import math
import numpy as np

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

class Molecule:
  def __init__(self, location, radius):
    self.location = location
    self.binding_partner = 0
    self.radius = radius
    self.locations = [self.location.copy()]
    self.R = get_R(np.random.uniform()*2*np.pi, np.random.uniform()*np.pi, np.random.uniform()*2*np.pi)
    self.dock_offset = np.matmul(np.array([0, 0, radius]), self.R) 
    self.dock_locations = [self.dock_offset + self.location]
  
  def move(self):
    D = 1.36e-6 / (100**2) * (10**10)**2# angstrom2/s
    dt = 10*1e-9 # s
    S = math.sqrt(2*D*dt)
    self.location[0] += S*np.random.normal()
    self.location[1] += S*np.random.normal()
    self.location[2] += S*np.random.normal()
    if np.any((self.location > 160)|(self.location < 0)):
      #self.location = np.random.randint(0, 160, 3).astype(float)
      self.location = np.array([50, 50, 50], float)
      self.dock_location = self.dock_offset + self.location
    self.locations.append(self.location.copy())

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
    self.dock_location = self.dock_offset + self.location
    self.dock_locations.append(self.dock_location)
  
  def set_id(self, id):
    self.id = id

  def plot(self, ax, t_step):
    limits = 100
    ax.set_xlim3d(-limits, limits)
    ax.set_ylim3d(-limits, limits)
    ax.set_zlim3d(-limits, limits)
    x = np.array([self.locations[t_step][0], self.dock_locations[t_step][0]])
    y = np.array([self.locations[t_step][1], self.dock_locations[t_step][1]])
    z = np.array([self.locations[t_step][2], self.dock_locations[t_step][2]])
    ax.plot3D(x, y, z)

    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v) * self.radius
    y = np.sin(u)*np.sin(v) * self.radius
    z = np.cos(v) * self.radius
    ax.plot_surface(x + self.locations[t_step][0], y + self.locations[t_step][1], z + self.locations[t_step][2], color='r', alpha=0.5)



class Ligand(Molecule):
  def hellow(self):
    print('hi')

class Substrate(Molecule):
  def hellow(self):
    print('hi')
  
  def calculate_ligand_distances(self, ligands):
    for ligand in ligands:
      distance = np.linalg.norm(self.dock_location - ligand.dock_location)
      if distance < 2:
        self.binding_partner = ligand.id
        ligand.binding_partner = self.id
      elif distance > 200:
        self.binding_partner = -1
        #ligand.location = np.array([50, 50, 50], float)
        #self.location = np.array([50, 50, 50+42], float)

        
