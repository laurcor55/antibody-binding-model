import math
import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import molecules as mol
import os
class Reaction:
  def __init__(self, molecules, rend):
    self.molecules = molecules
    self.end_reaction = False
    self.rend = rend
    self.dt = 1e-10
    self.binding = False
    id = 1
    for molecule in self.molecules:
      molecule.set_id(id)
      id += 1
    self.t_step = 1
    self.time = 0
    self.progress_reaction()

  def progress_reaction(self):
    while self.end_reaction == False:
      self.check_reaction_status()
      self.step_reaction()
      self.t_step += 1
      self.time += self.dt

  def check_reaction_status(self):
    substrates = [molecule for molecule in self.molecules if type(molecule) == mol.Substrate]
    ligands = [molecule for molecule in self.molecules if type(molecule) == mol.Ligand]
    for substrate in substrates:
      self.dt = substrate.calculate_ligand_distances(ligands, self.rend)
      if substrate.binding_partner > 0:
        self.end_reaction = True
        self.binding = True
        break
      elif substrate.binding_partner < 0:
        self.end_reaction = True
        self.binding = False
        break
 
  def show_animation(self):
    self.fig = plt.figure()
    self.ax = p3.Axes3D(self.fig)
    for t_step in range(self.t_step):
      self.ax.clear()
      for molecule in self.molecules:
        molecule.plot(self.ax, t_step)
      plt.pause(0.01)
    plt.show()

  def step_reaction(self):
    self.move_locked_molecules()
    self.move_unlocked_molecules()
    self.rotate_unlocked_molecules()
  
  def move_unlocked_molecules(self):
    overlap_check = False
    while overlap_check == False:
      for molecule in self.molecules:
        if molecule.locked_partner < 1:
          molecule.attempt_move(self.dt)
      overlap_check = self.check_overlaps()
    for molecule in self.molecules:
      if molecule.locked_partner < 1:
        molecule.move()
  
    
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

  def check_locks(self):
    ligands = [molecule for molecule in self.molecules if type(molecule) == mol.Ligand]
    for ligand in ligands:
      if ligand.locked_partner > 0:
        substrate = self.molecules[ligand.locked_partner-1]
        new_distances = [mol.calculate_distance(substrate_dock_location, ligand_dock_location) for substrate_dock_location, ligand_dock_location in zip(substrate.new_dock_locations, ligand.new_dock_locations)]
        close_docks = [new_distance < 3 for new_distance in new_distances]
        print(sum(close_docks))
        if sum(close_docks) < 2:
          return False
    return True

  def move_locked_molecules(self):
    ligands = [molecule for molecule in self.molecules if type(molecule) == mol.Ligand]
    for ligand in ligands:
      if ligand.locked_partner > 0:
        probability = np.exp(-4.2/(273+25)*ligand.radius*1e-10)
        substrate = self.molecules[ligand.locked_partner-1]
        number = np.random.uniform()
        if number < probability:
          ligand.location = find_locked_location(substrate, ligand)
          ligand.dock_locations = substrate.dock_locations
          ligand.R = ligand.location - ligand.dock_locations 
          
          ligand.dock_locations_over_time.append(ligand.dock_locations.copy())
          ligand.location_over_time.append(ligand.location.copy())
          substrate.dock_locations_over_time.append(substrate.dock_locations.copy())
          substrate.location_over_time.append(substrate.location.copy())
        else: 
          ligand.locked_partner = 0
          substrate.locked_partner = 0



def find_locked_location(substrate, ligand):
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
  new_location = np.array([x0, y0, z0]) - d * u
  return new_location

def scientific(input):
  output = "{:e}".format(input)
  return output
