import math
import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import molecules as mol
import os
from matplotlib.widgets import Slider
import time
class Reaction:
  def __init__(self, molecules, rend, minimum_binding_docks):
    self.reaction_status_history_goal = ['binding', 'free']
    self.molecules = molecules
    self.rend = rend
    self.minimum_binding_docks = minimum_binding_docks
    self.end_reaction = False
    self.binding_distance = 2
    self.time = 0
    self.calculate_dt()
    self.assign_molecule_id()
    self.check_molecule_status()
    self.reaction_status_history = [self.reaction_status]

  def assign_molecule_id(self):
    id = 1
    for molecule in self.molecules:
      molecule.set_id(id)
      id += 1

  def progress_reaction(self):
  #  while (self.end_reaction == False):
    while len(self.molecules[0].location_over_time)<10000:
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
    minimum_substrate_ligand_distance = self.rend*2
    for substrate in substrates:
      for ligand in ligands:
        distances = calculate_distance_docks(substrate.dock_locations, ligand.dock_locations)
        distance_min = min(distances)
        substrate.binding_partner = 0
        ligand.binding_partner = 0
        substrate.locked_partner = 0
        ligand.locked_partner = 0
        if distance_min < minimum_substrate_ligand_distance:
          self.minimum_substrate_ligand_distance = distance_min
        close_docks_num = sum(distance < self.binding_distance for distance in distances)
        if (close_docks_num >= self.minimum_binding_docks):
          substrate.binding_partner = ligand.id
          ligand.binding_partner = substrate.id
          self.reaction_status = 'binding'
        elif close_docks_num == 2:
          substrate.locked_partner = ligand.id
          ligand.locked_partner = substrate.id
        elif close_docks_num < 2:
          substrate.locked_partner = ligand.id
          ligand.locked_partner = substrate.id
          self.reaction_status = 'free'
    
  def check_reaction_status(self):
    if self.reaction_status_history[-1] != self.reaction_status:
      self.reaction_status_history.append(self.reaction_status)
    if self.reaction_status_history == self.reaction_status_history_goal:
      self.end_reaction = True
    elif self.minimum_substrate_ligand_distance > self.rend:
      self.end_reaction = True

  def calculate_dt(self):
    substrates = [molecule for molecule in self.molecules if (type(molecule) == mol.Substrate) or (type(molecule) == mol.FixedSubstrate)]
    ligands = [molecule for molecule in self.molecules if type(molecule) == mol.Ligand]
    distance_min_overall = self.rend*2
    for substrate in substrates:
      for ligand in ligands:
        distances = calculate_distance_docks(substrate.dock_locations, ligand.dock_locations)
        distance_min = min(distances)
        if distance_min < distance_min_overall:
          distance_min_overall = distance_min
          D_ligand = ligand.D
          D_substrate = substrate.D
    D = D_ligand + D_substrate

    distance = (distance_min_overall+1) *1e-10 # meters
    lamda = 1
    self.dt = 1/(12*D)*(distance/lamda)**2
 
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
    overlap_check = False
    while overlap_check == False:
      for molecule in self.molecules:
        if molecule.locked_partner < 1 and molecule.binding_partner < 1:
          molecule.attempt_move(self.dt)
        elif type(molecule) == mol.Ligand:
          if molecule.binding_partner > 0 :
            ligand = molecule
            substrate = self.molecules[molecule.binding_partner-1]
            self.find_bound_location(substrate, ligand)
          elif molecule.locked_partner > 0:
            ligand = molecule
            substrate = self.molecules[molecule.locked_partner-1]
            self.find_locked_location(substrate, ligand)
      overlap_check = self.check_overlaps()

  def check_overlaps(self):
    n = len(self.molecules)
    for ii in range(n-1):
      for jj in range(ii+1, n):
        if ii != jj:
          double_radius = self.molecules[ii].radius + self.molecules[jj].radius
          distance = np.linalg.norm(self.molecules[ii].new_location - self.molecules[jj].new_location)
          if distance < double_radius:
            return False
    return True

  def find_locked_location(self, substrate, ligand):
    new_location_check = False
    while new_location_check == False:
      substrate.attempt_move(self.dt)
      ligand.attempt_move(self.dt)
      distances = calculate_distance_docks(substrate.new_dock_locations, ligand.new_dock_locations)
      bound = np.sum(np.array(distances, float) <= self.binding_distance)
      if bound < 2:
        Poff = np.exp(-7.1)
        number = np.random.uniform()
        if number < Poff:
          new_location_check = True
      else:
        new_location_check = True
      
  def find_bound_location(self, substrate, ligand):
    new_location_check = False
    while new_location_check == False:
      substrate.attempt_move(self.dt)
      ligand.attempt_move(self.dt)
      distances = calculate_distance_docks(substrate.new_dock_locations, ligand.new_dock_locations)
      bound = np.sum(np.array(distances, float) <= self.binding_distance)
      if bound < self.minimum_binding_docks:
        Poff = np.exp(-11.5)
        number = np.random.uniform()
        if number < Poff:
          new_location_check = True
      else:
        new_location_check = True

 

def calculate_distance(location_1, location_2):
  return np.linalg.norm(location_1 - location_2) 

def calculate_distance_docks(dock_locations_1, dock_locations_2):
  distances = [calculate_distance(substrate_dock_location, ligand_dock_location) for substrate_dock_location, ligand_dock_location in zip(dock_locations_1, dock_locations_2)]
  return distances

def scientific(input):
  output = "{:e}".format(input)
  return output

def get_U(substrate_dock_locations, ligand_dock_locations, binding_distance):
  EI = 4.7e-24 # Nm^2
  E = EI/1.8e-32 #N/m^2
  klong = E*(2.7e-9*5.15e-9)/4e-9 #N/m
  klongi = klong
  kB=1.38e-23
  T=297
  distances = calculate_distance_docks(substrate_dock_locations, ligand_dock_locations)
  r = np.multiply(distances, 1e-10) 
  bound = np.array(distances, float) <= binding_distance
  return (1/2*klongi*r**2)/(kB*T)*bound 
