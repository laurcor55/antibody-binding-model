import math
import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import molecules as mol
import os
from matplotlib.widgets import Slider
import time
import numba as nb
from numba import jit


class Reaction:
  def __init__(self, ligand, substrates, rend, minimum_binding_docks):
    self.ligand = ligand
    self.substrates = substrates
    self.rend = rend
    self.minimum_binding_docks = minimum_binding_docks
    self.end_reaction = False
    self.binding_distance = 2
    self.time = 0
    self.reaction_status = 'free'
    self.check_molecule_status()
    self.calculate_dt()
  
  def back_to_start(self):
    self.end_reaction = False
    self.time = 0
    self.reaction_status = 'free'
    self.calculate_dt()
    [substrate.back_to_start() for substrate in self.substrates]
    self.ligand.back_to_start()
    self.check_molecule_status()

  def progress_reaction(self):
    while (self.end_reaction == False):
      self.calculate_dt()
      self.step_reaction()
      self.check_molecule_status()
      self.check_reaction_status()

  def check_molecule_status(self):
    self.minimum_substrate_ligand_dock_distance = self.rend*2
    for substrate in self.substrates:
      distances = calculate_distance_docks(substrate.dock_locations, self.ligand.dock_locations)
      distance_min = min(distances)
      if distance_min < self.minimum_substrate_ligand_dock_distance:
        self.minimum_substrate_ligand_dock_distance = distance_min
      close_docks_num = sum(distance < self.binding_distance for distance in distances)
      if (close_docks_num >= self.minimum_binding_docks):
        self.reaction_status = 'binding'
      if close_docks_num > 1:
        substrate.locked = True
      else:
        substrate.locked = False
  
  def check_reaction_status(self):
    if self.reaction_status == 'binding' or self.minimum_substrate_ligand_dock_distance > self.rend:
      self.end_reaction = True

  def calculate_dt(self):
    D_ligand = self.ligand.D
    D_substrate = self.substrates[0].D
    D = D_ligand + D_substrate
    distance = self.minimum_substrate_ligand_dock_distance *1e-10 # meters
    lamda = 5
    self.dt = 1/(12*D)*(distance/lamda)**2
 
  def show_animation(self):
    self.fig = plt.figure()
    self.slider_ax = plt.axes([0.15, 0.05, 0.75, 0.05])
    
    self.ax = p3.Axes3D(self.fig, [0, 0.1, 1, 1])
    t_steps = len(self.ligand.location_over_time)
    self.slider = Slider(self.slider_ax, label='t-step', valmin=0, valmax=t_steps-1, valinit=t_steps-1)
    self.slider.on_changed(self.update_slider)
    self.show_molecules_at_time(t_steps-1)
    plt.show()
  
  def update_slider(self, val):
    t_step = round(self.slider.val)
    self.show_molecules_at_time(t_step)

  def show_molecules_at_time(self, t_step):
    self.ax.clear()
    for substrate in self.substrates:
      substrate.plot(self.ax, t_step)
    self.ligand.plot(self.ax, t_step)
  
  def step_reaction(self):
    self.find_new_molecule_locations()
    self.move_molecules()
    self.time += self.dt
  
  def move_molecules(self):
    for substrate in self.substrates:
      substrate.move()
    self.ligand.move()
  
  def find_new_molecule_locations(self):
    no_overlaps_found = False
    while no_overlaps_found == False:
      ligand_locked = False
      for substrate in self.substrates:
        if substrate.locked == False:
          substrate.attempt_move(self.dt)
        else: 
          self.find_locked_location(substrate, self.ligand)
          ligand_locked = True
      if ligand_locked == False:
        self.ligand.attempt_move(self.dt)
      no_overlaps_found = self.check_overlaps()
    #if no_overlaps_found == False:
    #  for substrate in self.substrates:
    #    substrate.move_back()
    #  self.ligand.move_back()


  def check_overlaps(self):
    no_overlaps_found = True
    n = len(self.substrates)
    for ii in range(n-1):
      for jj in range(ii+1, n):
        if ii != jj:
          double_radius = self.substrates[ii].radius + self.substrates[jj].radius
          distance = np.linalg.norm(self.substrates[ii].new_location - self.substrates[jj].new_location)
          if distance < double_radius:
            no_overlaps_found = False
    for substrate in self.substrates:
      double_radius = substrate.radius + self.ligand.radius
      distance = np.linalg.norm(substrate.new_location - self.ligand.new_location)
      if distance < double_radius:
        no_overlaps_found = False
    return no_overlaps_found
  

  def find_locked_location(self, substrate, ligand):
    substrate.attempt_move(self.dt)
    ligand.attempt_move(self.dt)
    
    distances = calculate_distance_docks(substrate.new_dock_locations, ligand.new_dock_locations)
    new_bound_count = count_bound_docks(distances, self.binding_distance)

    if new_bound_count < 2:
      return_to_old_location = determine_thermodynamics(-7.1)
      if return_to_old_location:
        substrate.move_back()
        ligand.move_back()

  def __dict__(self):
    result = {'type': str(type(self)), 'rend': self.rend, 'reaction_status': self.reaction_status}
    return result


@jit(nopython=True)
def determine_thermodynamics(delta_G):
  Poff = np.exp(delta_G)
  number = np.random.uniform(0, 1, 1)
  if number < Poff:
    return_to_old_location = False
  else:
    return_to_old_location = True
  return return_to_old_location

@jit(nopython=True)
def calculate_distance_docks(dock_locations_1, dock_locations_2):
  distances = np.zeros(len(dock_locations_1))
  for ii in range(len(dock_locations_1)):
    distances[ii] = np.linalg.norm(dock_locations_1[ii, :]-dock_locations_2[ii, :])
  return distances

def scientific(input):
  output = "{:e}".format(input)
  return output

@jit(nopython=True)
def count_bound_docks(distances, binding_distance):
  bound = 0
  for distance in distances:
    if distance <= binding_distance:
      bound += 1
  return bound
