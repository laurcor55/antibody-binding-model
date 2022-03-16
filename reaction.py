import math
import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import molecules as mol

class Reaction:
  def __init__(self, molecules):
    self.molecules = molecules
    id = 1
    self.end_reaction = False
    for molecule in self.molecules:
      molecule.set_id(id)
      id += 1
    self.t_step = 0
    self.time = 0
    self.progress_reaction()
    self.show_animation()

  def progress_reaction(self):
    while self.end_reaction == False:
      self.step_reaction()
      self.compute_distances()
      self.t_step += 1
    print(self.t_step)

  def compute_distances(self):
    substrates = [molecule for molecule in self.molecules if type(molecule) == mol.Substrate]
    ligands = [molecule for molecule in self.molecules if type(molecule) == mol.Ligand]
    for substrate in substrates:
      substrate.calculate_ligand_distances(ligands)
      if substrate.binding_partner > 0:
        self.end_reaction = True
        break

  def show_animation(self):
    self.fig = plt.figure()
    self.ax = p3.Axes3D(self.fig)
    limits = 20
    if self.t_step < 500:
      for t_step in range(self.t_step):
        self.ax.clear()
        self.ax.set_xlim3d(-limits, limits)
        self.ax.set_ylim3d(-limits, limits)
        self.ax.set_zlim3d(-limits, limits)
        for molecule in self.molecules:
          molecule.plot(self.ax, t_step)
        plt.pause(0.01)
    else:
      for molecule in self.molecules:
        molecule.plot(self.ax, self.t_step)
    plt.show()
  
  def step_reaction(self):
    for molecule in self.molecules:
      molecule.move()
      molecule.rotate()

molecules = [mol.Ligand(np.array([0, 10, 10])), mol.Substrate(np.array([0, 0, 10])), mol.Ligand(np.array([0, 0, 0])), mol.Substrate(np.array([10, 0, 0]))]
Reaction(molecules)