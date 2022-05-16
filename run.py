import reaction as reac
import molecules as mol
import numpy as np
import time

radius = 18 # angstroms
minimum_binding_docks = 3
total_count = 1
rstart = 42 #angstroms
rend = 200 #angstroms
binding_count = 0
seed = int(time.time())
np.random.seed(seed)
times = []
time.perf_counter()
for kk in range(total_count):
#kk = 0
#while binding_count < 1:
  #molecules = [mol.Ligand([50, 50, 50], radius, [0, 0, 0]), mol.Substrate([50, 50, 50+rstart], radius, [np.random.uniform()*np.pi*2, np.random.uniform()*np.pi*2, np.random.uniform()*np.pi*2])]
  molecules = [mol.Ligand([50, 50, 50], radius, [0, 0, 0]), mol.Substrate([50, 50, 50+2*radius+1], radius, [0, np.pi, 0])]
 
  reaction = reac.Reaction(molecules, rend, minimum_binding_docks)
  reaction.progress_reaction()
  #binding_count += reaction.binding
  times.append(reaction.time)
  mean_kd = 1/np.mean(np.array(times))
  print(mean_kd)
  if kk % 10 ==0:
    print(str(kk) + ' of ' + str(total_count) + ', ' + str(binding_count) + ' bound', end="\r")
 # kk += 1
  del molecules
  del reaction

time.perf_counter()

#reaction.show_animation()
print(binding_count)
p = binding_count/total_count

print(p)
if False:
  D = molecules[0].D + molecules[1].D
  b = rstart*1e-10 # meters
  c = rend * 1e-10 #meters
  NA = 6.022e23

  k_d_b = 4*np.pi*b*D
  k_d_c = 4*np.pi*c*D
  k = k_d_b * p /(1-(1-p)*k_d_b/k_d_c) * 1000 * NA
  print(reac.scientific(k))
