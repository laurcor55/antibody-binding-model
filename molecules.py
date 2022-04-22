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
  
  return np.dot(np.dot(R_eta, R_phi), R_theta)

class Molecule:
  def __init__(self, location, radius, rotation):
    self.radius = radius
    self.location = np.array(location, float)
    
    temperature = 298
    boltzmann = 1.380649e-23
    viscocity = 8.9e-4
    self.D = boltzmann * temperature /(6*np.pi*viscocity*self.radius*1e-10) #m2/s
    self.D_r = boltzmann*temperature/(8*np.pi*viscocity*(self.radius*1e-10)**3)

    self.binding_partner = 0
    self.locked_partner = 0
    self.R = get_R(rotation[0], rotation[1], rotation[2])
    self.new_location = self.location
    self.location_over_time = [self.location.copy()]
    self.dock_offsets = [np.dot(dock_offset, self.R) for dock_offset in self.dock_offsets]
    self.dock_locations = [dock_offset + self.location for dock_offset in self.dock_offsets]
    self.dock_locations_over_time = [self.dock_locations.copy()]
  
  def __repr__(self):
    x = '__________\n' 
    x += 'center: ' + str(self.location) + '\n'
    x += 'dock locations: ' + str(self.dock_locations) + '\n'
    return x

  def move(self):
    self.location = self.new_location
    self.location_over_time.append(self.location.copy())    
  
  def set_id(self, id):
    self.id = id

  def plot(self, ax, t_step):
    limits = 100
    ax.set_xlim3d(-limits, limits)
    ax.set_ylim3d(-limits, limits)
    ax.set_zlim3d(-limits, limits)
    colors = ['b', 'g', 'r', 'y']
    ii = 0
    for dock_locations_over_time in self.dock_locations_over_time[t_step]:
      x = np.array([self.location_over_time[t_step][0], dock_locations_over_time[0]])
      y = np.array([self.location_over_time[t_step][1], dock_locations_over_time[1]])
      z = np.array([self.location_over_time[t_step][2], dock_locations_over_time[2]])
      ax.plot3D(x, y, z, color=colors[ii])
      ii += 1

    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v) * self.radius
    y = np.sin(u)*np.sin(v) * self.radius
    z = np.cos(v) * self.radius
    ax.plot_surface(x + self.location_over_time[t_step][0], y + self.location_over_time[t_step][1], z + self.location_over_time[t_step][2], color='r', alpha=0.5)

class Ligand(Molecule):
  def __init__(self, location, radius, rotation):
    self.dock_offsets = [np.array([-8.5, 8.5, radius]), np.array([8.5, 8.5, radius]), np.array([-8.5, -8.5, radius]), np.array([8.5, -8.5, radius])]
    super().__init__(location, radius, rotation)

  def attempt_move(self, dt):
    D = self.D * (10**10)**2# angstrom2/s
    S = math.sqrt(2*D*dt)
    nums = np.random.normal(scale=1, size=3)
    offset = np.multiply(S, nums)
    self.new_location = np.add(self.location, offset)
  
  def rotate(self, dt):
    S = math.sqrt(2*self.D_r*dt)
    delta_eta = S * np.random.normal()
    delta_phi = S * np.random.normal()
    delta_theta = S * np.random.normal()

    R = get_R(delta_theta, delta_phi, delta_eta)
    self.R = np.dot(self.R, R)
    self.dock_offsets = [np.dot(self.R, dock_offset) for dock_offset in self.dock_offsets]
    self.dock_locations = [dock_offset + self.location for dock_offset in self.dock_offsets]
    self.dock_locations_over_time.append(self.dock_locations.copy())
 

class Substrate(Molecule):
  def __init__(self, location, radius, rotation):
    self.dock_offsets =[np.array([8.5, 8.5, radius]), np.array([-8.5, 8.5, radius]), np.array([8.5, -8.5, radius]), np.array([-8.5, -8.5, radius])]
    super().__init__(location, radius, rotation)

  def attempt_move(self, dt):
    D = self.D * (10**10)**2# angstrom2/s
    S = math.sqrt(2*D*dt)
    nums = np.random.normal(scale=1, size=3)
    offset = np.multiply(S, nums)
    self.new_location = np.add(self.location, offset)
  
  def rotate(self, dt):
    S = math.sqrt(2*self.D_r*dt)
    delta_eta = S * np.random.normal()
    delta_phi = S * np.random.normal()
    delta_theta = S * np.random.normal()

    R = get_R(delta_theta, delta_phi, delta_eta)
    self.R = np.dot(self.R, R)
    self.dock_offsets = [np.dot(self.R, dock_offset) for dock_offset in self.dock_offsets]
    self.dock_locations = [dock_offset + self.location for dock_offset in self.dock_offsets]
    self.dock_locations_over_time.append(self.dock_locations.copy())


class FixedSubstrate(Molecule):
  def __init__(self, location, radius, rotation):
    self.dock_offsets =[np.array([8.5, 8.5, radius]), np.array([-8.5, 8.5, radius]), np.array([8.5, -8.5, radius]), np.array([-8.5, -8.5, radius])]
    super().__init__(location, radius, rotation)
 
  def attempt_move(self, dt):
    self.new_location = self.location
  
  def rotate(self, dt):
    self.dock_locations_over_time.append(self.dock_locations.copy())


def calculate_distance(location_1, location_2):
  return np.linalg.norm(location_1 - location_2) 
