import reaction as reac
import molecules as mol
import numpy as np
import time

radius = 18 # angstroms
minimum_binding_docks = 3
total_count = 10000
rstart = 42 #angstroms
rend = 200 #angstroms
binding_count = 0
seed = int(time.time())
np.random.seed(seed)
print(seed)
for kk in range(total_count):
#kk = 0
#while binding_count < 1:
  molecules = [mol.Ligand([50, 50, 50], radius, [0, 0, 0]), mol.Substrate([50, 50, 50+rstart], radius, [np.random.uniform()*np.pi*2, np.random.uniform()*np.pi*2, np.random.uniform()*np.pi*2])]
  #molecules = [mol.Substrate([50, 50, 50], radius, [0, 0, 0]), mol.Ligand([50, 78, 78], radius, [0, np.pi, -np.pi/2])]
  #molecules = [mol.Substrate([50, 50, 50], radius, [0, 0, 0]), mol.Ligand([50, 43, 87], radius, [0, np.pi, np.pi/8])]
  reaction = reac.Reaction(molecules, rend, minimum_binding_docks)
 # distances = reac.calculate_distance_docks(reaction.molecules[0].dock_locations, reaction.molecules[1].dock_locations)
  #print(distances)

 # print(reaction.molecules[0].locked_partner)
  reaction.progress_reaction()
  binding_count += reaction.binding
  if kk % 10 ==0:
    print(str(kk) + ' of ' + str(total_count) + ', ' + str(binding_count) + ' bound', end="\r")
 # kk += 1
  #del reaction
  #del molecules

reaction.show_animation()
print(binding_count)
p = binding_count/total_count

print(p)

D = molecules[0].D + molecules[1].D
b = rstart*1e-10 # meters
c = rend * 1e-10 #meters
NA = 6.022e23

k_d_b = 4*np.pi*b*D
k_d_c = 4*np.pi*c*D
k = k_d_b * p /(1-(1-p)*k_d_b/k_d_c) * 1000 * NA
print(reac.scientific(k))
