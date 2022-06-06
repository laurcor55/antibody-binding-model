#kernprof -l run.py
#lineprofilergui -l run.py.lprof
import reaction as reac
import molecules as mol
import numpy as np
import time

def calculate_kon(binding_count, total_count, molecules, rstart, rend):
  p = binding_count/total_count
  D = molecules[0].D + molecules[1].D
  b = rstart*1e-10 # meters
  c = rend * 1e-10 #meters
  NA = 6.022e23

  k_d_b = 4*np.pi*b*D
  k_d_c = 4*np.pi*c*D
  k = k_d_b * p /(1-(1-p)*k_d_b/k_d_c) * 1000 * NA
  return reac.scientific(k)

radius = 18 # angstroms
minimum_binding_docks = 3
total_count = 10000
rend = 200 #angstroms
binding_count = 0
seed = int(time.time())
np.random.seed(seed)

ligand_location = [0, 0, 39]
substrate_location = [0, 0, 0]
ligand_orientation = [0, np.pi, 0]
substrate_orientation = [0, 0, 0]
n_docks = 8

kk = 0
#while binding_count < 1:
while kk < total_count:
  molecules = [mol.Ligand(ligand_location, radius, ligand_orientation, n_docks), mol.FixedSubstrate(substrate_location, radius, substrate_orientation, n_docks)]
  rstart = np.linalg.norm(np.asarray(ligand_location) - np.asarray(substrate_location))
 
  reaction = reac.Reaction(molecules, rend, minimum_binding_docks, 'end_at_binding')
  reaction.progress_reaction()
  if reaction.reaction_status == 'binding':
    binding_count += 1
  if kk % 10 ==0:
    k_on = calculate_kon(binding_count, kk+1, molecules, rstart, rend)
    print(str(kk) + ' of ' + str(total_count) + ', ' + str(binding_count) + ' bound (' + reac.scientific(binding_count/(kk+1)) + '), ' +  k_on  + ' /M/s', end="\r")
  kk += 1
reaction.show_animation()
print(end='\n')
print(binding_count)
k_on = calculate_kon(binding_count, total_count, molecules)
print(k_on)