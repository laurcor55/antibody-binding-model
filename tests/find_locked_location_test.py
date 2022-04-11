import sys
sys.path.insert(1, '../')
import reaction
import molecules as mol
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import matplotlib.pyplot as plt

radius = 18 # angstroms
ligand = mol.Ligand(np.array([50, 50, 50], float), radius)
substrate = mol.Substrate(np.array([50+radius*2, 50, 50], float), radius)
ligand.new_location = reaction.find_locked_location(substrate, ligand)

ligand.dock_locations = substrate.dock_locations
ligand.R = ligand.new_location - ligand.dock_locations 
ligand.location = ligand.new_location

fig = plt.figure()
ax = p3.Axes3D(fig)


limits = 100
ax.set_xlim3d(-limits, limits)
ax.set_ylim3d(-limits, limits)
ax.set_zlim3d(-limits, limits)
colors = ['b', 'g', 'r', 'y']

ii = 0
for dock_location in ligand.dock_locations:
  x = np.array([ligand.location[0], dock_location[0]])
  y = np.array([ligand.location[1], dock_location[1]])
  z = np.array([ligand.location[2], dock_location[2]])
  ax.plot3D(x, y, z, color=colors[ii])
  ii+= 1

ii = 0
for dock_location in substrate.dock_locations:
  x = np.array([substrate.location[0], dock_location[0]])
  y = np.array([substrate.location[1], dock_location[1]])
  z = np.array([substrate.location[2], dock_location[2]])
  ax.plot3D(x, y, z, color=colors[ii])
  ii+= 1

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = np.cos(u)*np.sin(v) * substrate.radius
y = np.sin(u)*np.sin(v) * substrate.radius
z = np.cos(v) * substrate.radius
ax.plot_surface(x + substrate.location[0], y + substrate.location[1], z + substrate.location[2], color='r', alpha=0.5)

x = np.cos(u)*np.sin(v) * ligand.radius
y = np.sin(u)*np.sin(v) * ligand.radius
z = np.cos(v) * ligand.radius
ax.plot_surface(x + ligand.location[0], y + ligand.location[1], z + ligand.location[2], color='r', alpha=0.5)


plt.show()