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
  def __init__(self, molecules, rend, minimum_binding_docks, end_criteria):
    self.molecules = molecules
    self.rend = rend
    self.minimum_binding_docks = minimum_binding_docks
    self.end_reaction = False
    self.binding_distance = 2
    self.time = 0
    self.end_criteria = getattr(self, end_criteria)
    self.assign_molecule_id()
    self.reaction_status = 'free'
    self.calculate_dt()
    self.check_molecule_status()
    self.reaction_status_history = [self.reaction_status]

  def assign_molecule_id(self):
    id = 1
    for molecule in self.molecules:
      molecule.set_id(id)
      id += 1
  def progress_reaction(self):
    while (self.end_reaction == False):
   # while len(self.molecules[0].location_over_time)<1000:
      self.calculate_dt()
      self.step_reaction()
      self.check_molecule_status()
      self.check_reaction_status()

  def get_substrates(self):
    substrates = [molecule for molecule in self.molecules if (type(molecule) == mol.Substrate) or (type(molecule) == mol.FixedSubstrate)]
    return substrates

  def get_ligands(self):
    ligands = [molecule for molecule in self.molecules if type(molecule) == mol.Ligand]
    return ligands 

  def check_molecule_status(self):
    substrates = self.get_substrates()
    ligands = self.get_ligands()
    self.minimum_substrate_ligand_dock_distance = self.rend*2
    self.minimum_substrate_ligand_distance = self.rend*2
    for substrate in substrates:
      for ligand in ligands:
        distances = calculate_distance_docks_2(substrate.dock_locations, ligand.dock_locations)
        distance_min = min(distances)
        self.minimum_substrate_ligand_distance = np.linalg.norm(substrate.location - ligand.location)
        substrate.binding_partner = 0
        ligand.binding_partner = 0
        substrate.locked_partner = 0
        ligand.locked_partner = 0
        if distance_min < self.minimum_substrate_ligand_dock_distance:
          self.minimum_substrate_ligand_dock_distance = distance_min
        close_docks_num = sum(distance < self.binding_distance for distance in distances)
        if (close_docks_num >= self.minimum_binding_docks):
          substrate.binding_partner = ligand.id
          ligand.binding_partner = substrate.id
          self.reaction_status = 'binding'
        if close_docks_num == 0:
          self.reaction_status = 'free'
        if close_docks_num > 1:
          substrate.locked_partner = ligand.id
          ligand.locked_partner = substrate.id
    
  def end_at_distance(self):
    if self.minimum_substrate_ligand_dock_distance > self.rend:
      self.end_reaction = True
  
  def end_at_binding(self):
    if self.reaction_status == 'binding' or self.minimum_substrate_ligand_dock_distance > self.rend:
      self.end_reaction = True
  
  def end_at_unbinding(self):
    if self.reaction_status == 'free':
      self.end_reaction = True
  
    
  def check_reaction_status(self):
    if self.reaction_status_history[-1] != self.reaction_status:
      self.reaction_status_history.append(self.reaction_status)
    self.end_criteria()
    

  def calculate_dt(self):
    substrates = [molecule for molecule in self.molecules if (type(molecule) == mol.Substrate) or (type(molecule) == mol.FixedSubstrate)]
    ligands = [molecule for molecule in self.molecules if type(molecule) == mol.Ligand]
    distance_min_overall = self.rend*2
    D_ligand = ligands[0].D
    D_substrate = substrates[0].D
    for substrate in substrates:
      for ligand in ligands:
        distances = calculate_distance_docks_2(substrate.dock_locations, ligand.dock_locations)
        distance_min = min(distances)
        if distance_min < distance_min_overall:
          distance_min_overall = distance_min
          D_ligand = ligand.D
          D_substrate = substrate.D
    D = D_ligand + D_substrate

    distance = (distance_min_overall ) *1e-10 # meters
    lamda = 1
    self.dt = 1/(12*D)*(distance/lamda)**2
    if self.reaction_status == 'binding':
      self.dt = 1e-5
 
  def show_animation(self):
    self.fig = plt.figure()
    self.slider_ax = plt.axes([0.15, 0.05, 0.75, 0.05])
    
    self.ax = p3.Axes3D(self.fig, [0, 0.1, 1, 1])
    t_steps = len(self.molecules[0].location_over_time)
    self.slider = Slider(self.slider_ax, label='t-step', valmin=0, valmax=t_steps-1, valinit=t_steps-1)
    self.slider.on_changed(self.update_slider)
    self.show_molecules_at_time(t_steps-1)
    plt.show()
  
  def update_slider(self, val):
    t_step = round(self.slider.val)
    self.show_molecules_at_time(t_step)

  def show_molecules_at_time(self, t_step):
    self.ax.clear()
    for molecule in self.molecules:
      molecule.plot(self.ax, t_step)
  
  def step_reaction(self):
    self.find_new_molecule_locations()
    for molecule in self.molecules:
      molecule.move()
    self.time += self.dt
  
  def find_new_molecule_locations(self):
    no_overlaps_found = False
    while no_overlaps_found == False:
      for molecule in self.molecules:
        if molecule.locked_partner < 1 and molecule.binding_partner < 1:
          molecule.attempt_move(self.dt)
        elif molecule.binding_partner > 0:
          if type(molecule) == mol.Ligand:
            ligand = molecule
            substrate = self.molecules[molecule.binding_partner-1]
            self.find_bound_location_3(substrate, ligand)
        elif molecule.locked_partner > 0:
          if type(molecule) == mol.Ligand:
            ligand = molecule
            substrate = self.molecules[molecule.binding_partner-1]
            self.find_bound_location_2(substrate, ligand)
      no_overlaps_found = self.check_overlaps()

  def check_overlaps(self):
    no_overlaps_found = True
    n = len(self.molecules)
    for ii in range(n-1):
      for jj in range(ii+1, n):
        if ii != jj:
          double_radius = self.molecules[ii].radius + self.molecules[jj].radius
          distance = np.linalg.norm(self.molecules[ii].new_location - self.molecules[jj].new_location)
          if distance < double_radius:
            no_overlaps_found = False
    return no_overlaps_found
  
  
  def find_bound_location(self, substrate, ligand):
    substrate.attempt_move(self.dt)
    ligand.attempt_move(self.dt)

    U1 = get_U(substrate.dock_locations, ligand.dock_locations, self.binding_distance)
    U2 = get_U(substrate.new_dock_locations, ligand.new_dock_locations, self.binding_distance)
    dU = np.sum(U2 - U1)
    if dU > 0:
      delta_G = -1*dU
      return_to_old_location = determine_thermodynamics(delta_G, self.dt)
    else:
      return_to_old_location = False
    if return_to_old_location == True:
      substrate.move_back()
      ligand.move_back()

  def find_bound_location_2(self, substrate, ligand):
    distances = calculate_distance_docks_2(substrate.dock_locations, ligand.dock_locations)
    old_bound_count = count_bound_docks(distances, self.binding_distance)
    substrate.attempt_move(self.dt)
    ligand.attempt_move(self.dt)
    
    distances = calculate_distance_docks_2(substrate.new_dock_locations, ligand.new_dock_locations)
    new_bound_count = count_bound_docks(distances, self.binding_distance)

    if new_bound_count < old_bound_count:
      if old_bound_count == 2:
        delta_G = -7.1
      else:
        if new_bound_count == 2:
          delta_G = -11.5 
        else:
          delta_G = -18.6
      return_to_old_location = determine_thermodynamics(delta_G, self.dt)
      if return_to_old_location:
        substrate.move_back()
        ligand.move_back()
  
  def find_bound_location_3(self, substrate, ligand):
    distances = calculate_distance_docks_2(substrate.dock_locations, ligand.dock_locations)
    old_bound_count = count_bound_docks(distances, self.binding_distance)
    substrate.attempt_move(self.dt)
    ligand.attempt_move(self.dt)
    
    distances = calculate_distance_docks_2(substrate.new_dock_locations, ligand.new_dock_locations)
    new_bound_count = count_bound_docks(distances, self.binding_distance)

    if new_bound_count < old_bound_count:
      k_off = 0.2
      P_off = k_off*self.dt
      number = np.random.uniform(0, 1, 1)
      if number >= P_off:
        substrate.move_back()
        ligand.move_back()
        




@jit(nopython=True)
def determine_thermodynamics(delta_G, dt):
  Poff = np.exp(delta_G)
  number = np.random.uniform(0, 1, 1)
  if number < Poff:
    return_to_old_location = False
  else:
    return_to_old_location = True
  return return_to_old_location

@jit(nopython=True)
def calculate_distance_docks_2(dock_locations_1, dock_locations_2):
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

def get_U(substrate_dock_locations, ligand_dock_locations, binding_distance):
  EI = 4.7e-24 # Nm^2
  E = EI/1.8e-32 #N/m^2
  klong = E*(2.7e-9*5.15e-9)/4e-9 #N/m
  klongi = klong/4
 # klongi = 2e6*1e-7/4
  kB=1.38e-23
  T=297
  distances = calculate_distance_docks_2(substrate_dock_locations, ligand_dock_locations)
  r = np.multiply(distances, 1e-10) 
  bound = np.array(distances, float) <= binding_distance
  result = (klongi*r**2)/(kB*T)*bound 
  return result
