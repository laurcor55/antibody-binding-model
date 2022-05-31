import reaction as reac
import molecules as mol
import numpy as np
import time

end_criteria = lambda x:print(x)

radius = 18 # angstroms
minimum_binding_docks = 3
total_count = 10
rstart = 42 #angstroms
rend = radius*2+4 #angstroms
binding_count = 0
seed = int(time.time())
np.random.seed(seed)
times = []
unbinding_count = []


for kk in range(total_count):
  molecules = [mol.Ligand([50, 50, 50], radius, [0, 0, 0]), mol.Substrate([50, 50, 50+2*radius+1], radius, [0, np.pi, 0])]
 
  reaction = reac.Reaction(molecules, rend, minimum_binding_docks, 'end_at_distance')

  reaction.progress_reaction()
  
  unbinding_count.append(reaction.reaction_status_history.count('free'))
  times.append(reaction.time)
  mean_time = np.mean(np.array(times))
  mean_kd = 1/mean_time
  print(str(kk) + ' of ' + str(total_count) + ', ' + reac.scientific(mean_kd) + ' kd Mean unbinding events: ' + str(np.mean(unbinding_count)), end="\r")
  

print(end='\n')

reaction.show_animation()
