import math
import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import molecules as mol
import os
from matplotlib.widgets import Slider
class Reaction:
  def __init__(self, molecules, rend, minimum_binding_docks):
    self.molecules = molecules
    self.rend = rend
    self.minimum_binding_docks = minimum_binding_docks
    self.end_reaction = False
    self.binding = False
    self.binding_distance = 5
    self.assign_molecule_id()
    self.check_reaction_status()
    self.calculate_dt()

  def assign_molecule_id(self):
    id = 1
    for molecule in self.molecules:
      molecule.set_id(id)
      id += 1

  def progress_reaction(self):
    while self.end_reaction == False:
      self.calculate_dt()
      self.step_reaction()
      self.check_reaction_status()

  def get_substrates(self):
    substrates = [molecule for molecule in self.molecules if (type(molecule) == mol.Substrate) or (type(molecule) == mol.FixedSubstrate)]
    return substrates

  def get_ligands(self):
    ligands = [molecule for molecule in self.molecules if type(molecule) == mol.Ligand]
    return ligands 

  def check_reaction_status(self):
    substrates = self.get_substrates()
    ligands = self.get_ligands()
    minimum_substrate_ligand_distance = self.rend*2
    for substrate in substrates:
      for ligand in ligands:
        distances = [mol.calculate_distance(substrate_dock_location, ligand_dock_location) for substrate_dock_location, ligand_dock_location in zip(substrate.dock_locations, ligand.dock_locations)]
        distance_min = min(distances)
        if distance_min < minimum_substrate_ligand_distance:
          minimum_substrate_ligand_distance = distance_min
        close_docks_num = sum(distance < self.binding_distance for distance in distances)
        if close_docks_num >= self.minimum_binding_docks:
          substrate.binding_partner = ligand.id
          ligand.binding_partner = substrate.id
          self.end_reaction = True
          self.binding = True
        if close_docks_num >= 2:
          substrate.locked_partner = ligand.id
          ligand.locked_partner = substrate.id
    if minimum_substrate_ligand_distance > self.rend:
      self.end_reaction = True
      self.binding = False

  def calculate_dt(self):
    substrates = [molecule for molecule in self.molecules if (type(molecule) == mol.Substrate) or (type(molecule) == mol.FixedSubstrate)]
    ligands = [molecule for molecule in self.molecules if type(molecule) == mol.Ligand]
    distance_min_overall = self.rend*2
    for substrate in substrates:
      for ligand in ligands:
        distances = [calculate_distance(substrate_dock_location, ligand_dock_location) for substrate_dock_location, ligand_dock_location in zip(substrate.dock_locations, ligand.dock_locations)]
        distance_min = min(distances)
        if distance_min < distance_min_overall:
          distance_min_overall = distance_min
          D_ligand = ligand.D
          D_substrate = substrate.D
    D = D_ligand + D_substrate
    distance = distance_min_overall *1e-10 # meters
    lamda = 10 
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
    
  
  def move_unlocked_molecules(self):
    overlap_check = False
    while overlap_check == False:
      for molecule in self.molecules:
        if molecule.locked_partner < 1:
          molecule.attempt_move(self.dt)
      overlap_check = self.check_overlaps()

  def find_new_molecule_locations(self):
    overlap_check = False
    while overlap_check == False:
      for molecule in self.molecules:
        if molecule.locked_partner < 1:
          molecule.attempt_move(self.dt)
        elif type(molecule) == mol.Ligand:
          ligand = molecule
          substrate = self.molecules[molecule.locked_partner-1]
          self.find_locked_location_2(substrate, ligand)
      overlap_check = self.check_overlaps()

  def rotate_unlocked_molecules(self):
    for molecule in self.molecules:
      if molecule.locked_partner < 1:
        molecule.rotate(self.dt)

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

  def find_locked_location_2(self, substrate, ligand):
    new_location_check = False
    while new_location_check == False:
      substrate.attempt_move(self.dt)
      ligand.attempt_move(self.dt)
      U1 = get_U(substrate, ligand, self.binding_distance)
      U2 = get_new_U(substrate, ligand, self.binding_distance)
      dU = np.sum(U2 - U1)
      if dU > 0:
        Poff = np.exp(-dU)
        number = np.random.uniform()
        if number < Poff:
          new_location_check = True
      else:
        new_location_check = True
   #   probability = np.exp(-7.1)
   #   number = np.random.uniform()
   #   if number > probability:
   #     self.find_locked_location(substrate, ligand)
   #   else: 
   #     ligand.locked_partner = 0
   #     substrate.locked_partner = 0

  def unbind_molecules(self):
    bound_pairs = self.get_bound_pairs()
    for bound_pair in bound_pairs:
      substrate = bound_pair[0]
      ligand = bound_pair[1]
      probability = np.exp(-18.6)
      number = np.random.uniform()
      if number > probability:
        print('rotate?')
      else: 
        ligand.bound_partner = 0
        substrate.bound_partner = 0

  def get_locked_pairs(self):
    locked_pairs = []
    ligands = [molecule for molecule in self.molecules if type(molecule) == mol.Ligand]
    for ligand in ligands:
      if ligand.locked_partner > 0:
        substrate = self.molecules[ligand.locked_partner-1]
        locked_pairs.append([substrate, ligand])
    return locked_pairs
  
    def get_bound_pairs(self):
      bound_pairs = []
      ligands = [molecule for molecule in self.molecules if type(molecule) == mol.Ligand]
      for ligand in ligands:
        if ligand.bound_partner > 0:
          substrate = self.molecules[ligand.binding_partner-1]
          bound_pairs.append([substrate, ligand])
      return bound_pairs
  
  def find_locked_location(self, substrate, ligand):
    dock_locations = substrate.dock_locations
    p1 = dock_locations[0]
    p2 = dock_locations[1]
    p3 = dock_locations[2]
    p4 = dock_locations[3]

    x0 = (p1[0] + p2[0] +p3[0] + p4[0])/4
    y0 = (p1[1] + p2[1] +p3[1] + p4[1])/4
    z0 = (p1[2] + p2[2] +p3[2] + p4[2])/4

    x1 = substrate.location[0]
    y1 = substrate.location[1]
    z1 = substrate.location[2]
    
    v = np.array([x1, y1, z1]) - np.array([x0, y0, z0])
    u = v/np.linalg.norm(v)
    d = ligand.radius
    ligand.new_location = np.array([x0, y0, z0]) - d * u
    ligand.location = ligand.new_location
    ligand.dock_locations = substrate.dock_locations
    ligand.R = ligand.location - ligand.dock_locations 
    ligand.dock_locations_over_time.append(ligand.dock_locations.copy())
    ligand.location_over_time.append(ligand.location.copy())
    
    substrate.dock_locations_over_time.append(substrate.dock_locations.copy())
    substrate.location_over_time.append(substrate.location.copy())

def calculate_distance(location_1, location_2):
  return np.linalg.norm(location_1 - location_2) 

def scientific(input):
  output = "{:e}".format(input)
  return output

def get_U(substrate, ligand, binding_distance):
  EI = 4.7e-24 # Nm^2
  E = EI/1.8e-32 #N/m^2
  klong = E*(2.7e-9*5.15e-9)/4e-9 #N/m
  klongi = klong/3
  kB=1.38e-23
  T=297
  distances = [calculate_distance(substrate_dock_location, ligand_dock_location) for substrate_dock_location, ligand_dock_location in zip(substrate.dock_locations, ligand.dock_locations)]
  r = np.multiply(distances, 1e-10)
  bound = np.array(distances, float) <= binding_distance
  return (1/2*klongi*r**2)/(kB*T)*bound 

def get_new_U(substrate, ligand, binding_distance):
  EI = 4.7e-24 # Nm^2
  E = EI/1.8e-32 #N/m^2
  klong = E*(2.7e-9*5.15e-9)/4e-9 #N/m
  klongi = klong/3
  kB=1.38e-23
  T=297
  distances = [calculate_distance(substrate_dock_location, ligand_dock_location) for substrate_dock_location, ligand_dock_location in zip(substrate.new_dock_locations, ligand.new_dock_locations)]
  r = np.multiply(distances, 1e-10)
  bound = np.array(distances, float) <= binding_distance
  return (1/2*klongi*r**2)/(kB*T)*bound 