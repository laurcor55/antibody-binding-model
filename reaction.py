import math
import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import molecules as mol
import os

class Reaction:
  def __init__(self, molecules):
    self.molecules = molecules
    id = 1
    self.end_reaction = False
    self.binding = False
    for molecule in self.molecules:
      molecule.set_id(id)
      id += 1
    self.t_step = 0
    self.time = 0
    self.progress_reaction()
    #self.show_animation()

  def progress_reaction(self):
    while self.end_reaction == False:
   # while self.end_reaction == False and self.time < 0.000154734999999959:
      self.step_reaction()
      self.compute_distances()
      self.t_step += 1
      self.time += 10e-9

  def compute_distances(self):
    substrates = [molecule for molecule in self.molecules if type(molecule) == mol.Substrate]
    ligands = [molecule for molecule in self.molecules if type(molecule) == mol.Ligand]
    for substrate in substrates:
      substrate.calculate_ligand_distances(ligands)
      if substrate.binding_partner > 0:
        self.end_reaction = True
        self.binding = True
      elif substrate.binding_partner < 0:
        self.end_reaction = True
        self.binding = False
        break
 

  def show_animation(self):
    self.fig = plt.figure()
    self.ax = p3.Axes3D(self.fig)
    if self.t_step < 50000:
      for t_step in range(self.t_step):
        self.ax.clear()
        for molecule in self.molecules:
          molecule.plot(self.ax, t_step)
        plt.pause(0.1)
    else:
      for molecule in self.molecules:
        molecule.plot(self.ax, self.t_step)
    plt.show()
    plt.savefig('output.png')

  def step_reaction(self):
    for molecule in self.molecules:
      molecule.move()
      molecule.rotate()
    distance = np.linalg.norm(self.molecules[0].location - self.molecules[1].location)
    if distance < 36:
      for molecule in self.molecules:
        molecule.location = molecule.locations[-2]
        molecule.dock_location = molecule.dock_locations[-2]
      self.time += -10e-9
      self.t_step += -1

def scientific(input):
  output = "{:e}".format(input)
  return output

radius = 18
binding_count = 0
total_count = 500
rstart = 42
total_time = np.zeros(total_count)
binding = np.zeros(total_count)
for ii in range(total_count):
  molecules = [mol.Ligand(np.array([50, 50, 50], float), radius), mol.Substrate(np.array([50, 50, 50+rstart], float), radius)]
  reaction = Reaction(molecules)
  total_time[ii] = reaction.time
  binding[ii] = reaction.binding
  print(str(ii) + ' of ' + str(total_count))
p = np.sum(binding)/len(binding)
print(p)
print(binding)
print(total_time)
plt.hist(total_time)
plt.show()
mean_time = np.mean(total_time)
print(mean_time)
volume = 4/3*np.pi*(100e-10)**3 #(160e-10)**3 * 1000
moles = 1/6.022e23
concentration = moles/volume
#print(concentration)
#k = total_count/(concentration*total_time)
#print(scientific(k))
#print(scientific(1/(mean_time*concentration)))
D = 1.36e-6 / (100**2) 
b = 42e-10
c = 200e-10
k_d_b = 4*np.pi*b*D
k_d_c = 4*np.pi*c*D
NA = 6.022e23
k = k_d_b * p /(1-(1-p)*k_d_b/k_d_c) * 1000*NA
print(scientific(k))

#if p < 1:
#  lam = np.log(1-p)/(-0.000154734999999959)
#  t = 1/lam
#  k = 1/(t*concentration)
#  print(scientific(k))