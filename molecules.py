import math
import numpy as np
import scipy
from scipy.spatial.transform import Rotation as R
from numba import jit
import time



class Molecule:
  def __init__(self, location, radius, rotation):
    self.radius = radius
    self.location = np.array(location, float)
    
    temperature = 298
    boltzmann = 1.380649e-23
    viscocity = 8.9e-4
    self.D = boltzmann * temperature /(6*np.pi*viscocity*self.radius*1e-10) #m2/s
    self.D_r = boltzmann * temperature/(8*np.pi*viscocity*(self.radius*1e-10)**3)

    self.binding_partner = 0
    self.locked_partner = 0
    self.R = get_R(rotation[0], rotation[1], rotation[2])
    self.new_location = self.location
    self.location_over_time = [self.location]
    self.dock_offsets = multiply_R(self.R, self.dock_offsets)
    self.dock_locations = get_dock_locations(self.dock_offsets, self.location)
    self.dock_locations_over_time = [self.dock_locations]
  
  def __repr__(self):
    x = '__________\n' 
    x += 'center: ' + str(self.location) + '\n'
    x += 'dock locations: ' + str(self.dock_locations) + '\n'
    return x

  def attempt_move(self, dt):
    self.new_location = find_new_location(self.D, dt, self.location)
    delta = get_new_angles(self.D_r, dt)
    self.new_R = get_R(delta[0], delta[1], delta[2])
    self.new_dock_offsets = multiply_R(self.new_R, self.dock_offsets)
    self.new_dock_locations = get_dock_locations(self.new_dock_offsets, self.new_location)

  def move(self):
    self.location = self.new_location
    self.R = self.new_R
    self.dock_offsets = self.new_dock_offsets
    self.dock_locations = self.new_dock_locations
    self.location_over_time.append(self.location.copy())    
    self.dock_locations_over_time.append(self.dock_locations.copy())
  
  def move_back(self):
    self.new_location = self.location
    self.new_R = self.R
    self.new_dock_offsets = self.dock_offsets
    self.new_dock_locations = self.dock_locations

  
  def set_id(self, id):
    self.id = id

  def plot(self, ax, t_step):
    limits = 100
    ax.set_xlim3d(-limits, limits)
    ax.set_ylim3d(-limits, limits)
    ax.set_zlim3d(-limits, limits)
    colors = ['b', 'g', 'r', 'c', 'b', 'g', 'r', 'c']
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
    ax.plot_surface(x + self.location_over_time[t_step][0], y + self.location_over_time[t_step][1], z + self.location_over_time[t_step][2], color='r', alpha=0.3)

class Ligand(Molecule):
  def __init__(self, location, radius, rotation):
    dock_offsets_1 = np.array(([-8.5, 8.5, radius], [8.5, 8.5, radius], [-8.5, -8.5, radius], [8.5, -8.5, radius]))
    dock_offsets_2 = np.array(([8.5, -8.5, -1*radius], [8.5, 8.5, -1*radius], [-8.5, -8.5, -1*radius]))
    self.dock_offsets = dock_offsets_1
    #self.dock_offsets = np.concatenate((dock_offsets_1, dock_offsets_2))
    
    
    super().__init__(location, radius, rotation)

class Substrate(Molecule):
  def __init__(self, location, radius, rotation):
    dock_offsets_1 = np.array(([8.5, 8.5, radius], [-8.5, 8.5, radius], [8.5, -8.5, radius], [-8.5, -8.5, radius]))
    self.dock_offsets = dock_offsets_1
    #self.dock_offsets = np.concatenate((dock_offsets_1, dock_offsets_1))
    super().__init__(location, radius, rotation)

class FixedSubstrate(Molecule):
  def __init__(self, location, radius, rotation):
    self.dock_offsets = np.array(([8.5, 8.5, radius], [-8.5, 8.5, radius], [8.5, -8.5, radius]))#, np.array([-8.5, -8.5, radius])]
    super().__init__(location, radius, rotation)

  def attempt_move(self, dt):
    self.new_location = self.location
    self.new_R = self.R
    self.new_dock_offsets = self.dock_offsets
    self.new_dock_locations = self.dock_locations


def calculate_distance(location_1, location_2):
  return np.linalg.norm(location_1 - location_2) 

@jit(nopython=True)
def multiply_R(R, dock_offsets):
  nrows = dock_offsets.shape[0]
  new_dock_offsets = np.zeros((nrows, 3))
  for ii in range(nrows):
    dock_offsets_single = dock_offsets[ii, :]
    ra, ca = R.shape
    output = np.zeros(3)
    for i in range(ra):
      for k in range(3):
          output[i] += R[i, k] * dock_offsets_single[k]
    new_dock_offsets[ii, :] = output
  return new_dock_offsets

@jit(nopython=True)
def find_new_location(D, dt, location):
  S = np.sqrt(2*D*dt)*10**10
  new_location = np.zeros(3)
  for ii in range(3):
    new_location[ii] = np.random.normal() * S + location[ii]
  return new_location

@jit(nopython=True)
def get_new_angles(D_r, dt):
  S = np.sqrt(2*D_r*dt)
  delta = np.zeros(3)
  for ii in range(3):
    delta[ii] = S * np.random.normal()
  return delta

@jit(nopython=True)
def get_dock_locations(offsets, location):
  dock_locations = np.empty_like(offsets)
  for ii in range(dock_locations.shape[0]):
    dock_locations[ii] = offsets[ii] + location
  return dock_locations

@jit(nopython=True)
def get_R(theta, phi, eta):
  R_theta = np.array(((np.cos(theta), -np.sin(theta), 0.0), \
    (np.sin(theta), np.cos(theta), 0.0), \
    (0.0, 0.0, 1.0)))

  R_phi = np.array(((np.cos(phi), 0.0, np.sin(phi)), \
    (0.0, 1.0, 0.0), \
    (-np.sin(phi), 0.0, np.cos(phi))))

  R_eta = np.array(((1.0, 0.0, 0.0), \
    (0.0, np.cos(eta), -np.sin(eta)), \
    (0.0, np.sin(eta), np.cos(eta))))
  first_part = np.dot(R_eta, R_phi)
  return np.dot(first_part, R_theta)