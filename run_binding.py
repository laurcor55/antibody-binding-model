#kernprof -l run.py
#lineprofilergui -l run.py.lprof
import reaction as reac
import molecules as mol
import numpy as np
import time

radius = 18 # angstroms
minimum_binding_docks = 3
total_count = 100000
rstart = 42 #angstroms
rend = 200 #angstroms
binding_count = 0
seed = int(time.time())
np.random.seed(seed)

kk = 0
#while binding_count < 1:
while kk < total_count:
  molecules = [mol.Ligand([50, 50, 50], radius, np.random.uniform(size=3)*np.pi*2), mol.Substrate([50, 50, 50+rstart], radius, np.random.uniform(size=3)*np.pi*2)]
 
  reaction = reac.Reaction(molecules, rend, minimum_binding_docks, 'end_at_binding')

  reaction.progress_reaction()
  if reaction.reaction_status == 'binding':
    binding_count += 1
  if kk % 10 ==0:
    print(str(kk) + ' of ' + str(total_count) + ', ' + str(binding_count) + ' bound', end="\r")
  kk += 1
reaction.show_animation()
print(end='\n')
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
