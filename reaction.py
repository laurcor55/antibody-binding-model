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
    self.binding = 0
    id = 1
    for molecule in self.molecules:
      molecule.set_id(id)
      id += 1
    self.t_step = 0
    self.time = 0
    self.progress_reaction()

  def progress_reaction(self):
    self.dt = 1e-10
    while self.end_reaction == False:
      self.step_reaction()
      self.compute_distances()
      self.t_step += 1
      self.time += self.dt

  def compute_distances(self):
    substrates = [molecule for molecule in self.molecules if type(molecule) == mol.Substrate]
    ligands = [molecule for molecule in self.molecules if type(molecule) == mol.Ligand]
    for substrate in substrates:
      self.dt = substrate.calculate_ligand_distances(ligands, self.rend)
      if substrate.binding_partner > 0:
        self.end_reaction = True
        self.binding = 1
        break
      elif substrate.binding_partner < 0:
        self.end_reaction = True
        self.binding = 0
        break
 
  def show_animation(self):
    self.fig = plt.figure()
    self.ax = p3.Axes3D(self.fig)
    for t_step in range(self.t_step):
      self.ax.clear()
      for molecule in self.molecules:
        molecule.plot(self.ax, t_step)
      plt.pause(0.1)
    plt.show()

  def step_reaction(self):
    self.find_new_location()
    for molecule in self.molecules:
      molecule.move()
      molecule.rotate(self.dt)
    
  def find_new_location(self):
    self.overlap = True
    while self.overlap == True:
      for molecule in self.molecules:
        molecule.attempt_move(self.dt)
      self.check_overlaps()
    

  def check_overlaps(self):
    n = len(self.molecules)
    for ii in range(n-1):
      for jj in range(ii+1, n):
        if ii != jj:
          double_radius = self.molecules[ii].radius + self.molecules[jj].radius
          distance = np.linalg.norm(self.molecules[ii].new_location - self.molecules[jj].new_location)
          if distance < double_radius:
            self.overlap = True
            break
          #while distance < double_radius:
          #  self.molecules[ii].new_location += (self.molecules[ii].location - self.molecules[ii].#new_location)/10
          #  self.molecules[jj].new_location += (self.molecules[jj].location - self.molecules[jj].new_location)/10
          #  distance = np.linalg.norm(self.molecules[ii].new_location - self.molecules[jj].new_location)
    self.overlap = False
            

def scientific(input):
  output = "{:e}".format(input)
  return output

radius = 18
total_count = 5000
rstart = 42
rend = 200
binding_count = 0
time = np.zeros(total_count)
for kk in range(total_count):
  molecules = [mol.Ligand(np.array([50, 50, 50], float), radius), mol.Substrate(np.array([50, 50, 50+rstart], float), radius)]
  reaction = Reaction(molecules, rend)
  binding_count += reaction.binding
  time[kk] = reaction.time
  if kk % 10 ==0:
    print(str(kk) + ' of ' + str(total_count))
#reaction.show_animation()
print(binding_count)
p = binding_count/total_count
print(p)

D = 1.36e-6 / (100**2)
b = rstart*1e-10
c = rend * 1e-10
NA = 6.022e23

k_d_b = 4*np.pi*b*D
k_d_c = 4*np.pi*c*D
k = k_d_b * p /(1-(1-p)*k_d_b/k_d_c) * 1000 * NA
print(scientific(k))
