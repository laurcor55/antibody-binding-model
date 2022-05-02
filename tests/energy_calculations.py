import sys
sys.path.insert(1, '../')
import reaction as reac
import molecules as mol
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import matplotlib.pyplot as plt

def get_U(substrate, ligand):
  EI = 4.7e-24 # Nm^2
  E = EI/1.8e-32 #N/m^2
  klong = E*(2.7e-9*5.15e-9)/4e-9 #N/m
  klongi = klong/4
  print(klongi)
  kB=1.38e-23
  T=297
  distances = [reac.calculate_distance(substrate_dock_location, ligand_dock_location) for substrate_dock_location, ligand_dock_location in zip(substrate.dock_locations, ligand.dock_locations)]
  r = np.multiply(distances, 1e-10)
  bound = np.array(distances, float) <= 2
  return (1/2*klongi*r**2)/(kB*T)*bound 


molecules = [mol.Substrate([50, 50, 50], 9, [0, 0, 0]), mol.Ligand([50, 66, 66], 7, [0, np.pi, -np.pi/2])]


rend = 1000
minimum_binding_docks = 3
reaction = reac.Reaction(molecules, rend, minimum_binding_docks)
U1 = get_U(reaction.molecules[0], reaction.molecules[1])

reaction.molecules[0].location = np.array([50.1, 50.1, 50.1], float)

reaction.molecules[0].dock_locations = [np.add(0.1, dock_location) for dock_location in reaction.molecules[0].dock_locations]
U2 = get_U(reaction.molecules[0], reaction.molecules[1])
dU = np.sum(U2 - U1)
if dU > 0:
  Poff = np.exp(-dU)
  number = np.random.uniform()
  if number < Poff:
    print('allow')
#checkr=sqrt((alphax-PFbetax(1,:)).^2+(alphay-PFbetay(1,:)).^2+(alphaz-PFbetaz(1,:)).^2)
#checkrlat = sqrt(sum((Nlat-Dlat).^2))
#U2 = [(1/2*klongi*checkr.^2)/(kB*T).*bound,(1/2*klati*checkrlat.^2)/(kB*T).*latbound]

#dU = sum(U2-U1)