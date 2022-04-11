import reaction as reac
import molecules as mol
import numpy as np


radius = 18 # angstroms
total_count = 50000
rstart = 42 #angstroms
rend = 200 #angstroms
binding_count = 0
for kk in range(total_count):
#kk = 0
#while binding_count < 1:
  molecules = [mol.Ligand(np.array([50, 50, 50], float), radius), mol.Substrate(np.array([50, 50, 50+rstart], float), radius)]
  reaction = reac.Reaction(molecules, rend)
  binding_count += reaction.binding
  if kk % 10 ==0:
    print(str(kk) + ' of ' + str(total_count) + ', ' + str(binding_count) + ' bound')
  #kk += 1

#reaction.show_animation()
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
