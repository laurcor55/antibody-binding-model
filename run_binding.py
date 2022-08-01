#kernprof -l run.py
#lineprofilergui -l run.py.lprof
import reaction as reac
import molecules as mol
import numpy as np
import time
import json
import copy
import os
def append_output(new_data, file_name):
  with open(file_name, 'r+') as file:
    file_data = json.load(file)
    file_data.append(new_data)
    file.seek(0)
    json.dump(file_data, file, indent=4)

def write_output(output_dict, file_name):
  with open(file_name, "w") as outfile:
    json.dump([output_dict], outfile)

def to_dict(ligand_dict, substrates_dict, reaction_dict):
  overall_json = {}
  overall_json['reaction'] = reaction_dict
  overall_json['ligand'] = ligand_dict
  overall_json['substrates'] = substrates_dict

  return overall_json

def calculate_kon(binding_count, total_count, ligand, substrates, rstart, rend):
  p = binding_count/total_count
  D = substrates[0].D + ligand.D
  b = rstart*1e-10 # meters
  c = rend * 1e-10 #meters
  NA = 6.022e23
  k_d_b = 4*np.pi*b*D
  k_d_c = 4*np.pi*c*D
  k = k_d_b * p /(1-(1-p)*k_d_b/k_d_c) * 1000 * NA
  return reac.scientific(k)

def run_reactions(file_name, substrates, ligand):
  minimum_binding_docks = 3
  total_count = 1
  rend = 200 #angstroms
  binding_count = 0
  seed = int(time.time())
  np.random.seed(seed)
  substrates_dict = [molecule.status() for molecule in substrates]
  ligand_dict = ligand.status()
  rstart = np.linalg.norm(np.asarray(ligand.location) - np.asarray(substrates[0].location))
  kk = 0
  
  reaction = reac.Reaction(ligand, substrates, rend, minimum_binding_docks)
  center_binding_count = 0
  
  while binding_count < 40:
  #while kk < total_count:
    reaction.back_to_start()
    reaction.progress_reaction()
   # reaction.show_animation()
    if reaction.reaction_status == 'binding':
      binding_count += 1
      if (reaction.substrates[0].locked == True):
        center_binding_count += 1
    if kk % 10 ==0:
      k_on = calculate_kon(binding_count, kk+1, ligand, substrates, rstart, rend)
      print(str(kk) + ' of ' + str(total_count) + ', ' + str(binding_count) + ' bound (' + reac.scientific(binding_count/(kk+1)) + '), ' +  k_on  + ' /M/s', end="\r")
    kk += 1
  total_count = kk
  center_binding_proportion = center_binding_count/binding_count

  reaction_dict = {'total_count': total_count, 'binding_count': binding_count, 'minimum_binding_docks': minimum_binding_docks, 'center_binding_proportion': center_binding_proportion}
  output_list = to_dict(ligand_dict, substrates_dict, reaction_dict)
  
  print(end='\n')
  print(binding_count)
  k_on = calculate_kon(binding_count, total_count, ligand, substrates, rstart, rend)
  print(k_on)
  if os.path.exists(file_name):
    append_output(output_list, file_name)
  else:
    write_output(output_list, file_name)

if __name__ == "__main__":
  file_name = 'figures/7-8-22/output_all_random.json'

  radius = 18 # angstroms
  z_spacings = [42]
  degs = [0]
  substrate_spacing = 36
  dock_rotations = [[0, 0, 0]]#, [0, np.pi/4, 0]]

  for z_spacing in z_spacings:
    for deg in degs:
      ligand_location = [0, 0, z_spacing]
      substrate_location = [0, 0, 0]
      substrate_location_2 = [0, substrate_spacing, 0]
      substrate_location_3 = [0, -substrate_spacing, 0]
      substrate_location_4 = [substrate_spacing, 0, 0]
      substrate_location_5 = [-substrate_spacing, 0, 0]
      substrate_location_6 = [substrate_spacing, substrate_spacing, 0]
      substrate_location_7 = [-substrate_spacing, -substrate_spacing, 0]
      substrate_location_8 = [-substrate_spacing, substrate_spacing, 0]
      substrate_location_9 = [substrate_spacing, -substrate_spacing, 0]

      

      substrates_start = [mol.FixedSubstrate(substrate_location, radius, [np.random.rand()*np.pi*2, 0, 0], dock_rotations)]#, mol.FixedSubstrate(substrate_location_2, radius, [np.random.rand()*np.pi*2, 0, 0], dock_rotations), mol.FixedSubstrate(substrate_location_3, radius, [np.random.rand()*np.pi*2, 0, 0], dock_rotations), mol.FixedSubstrate(substrate_location_4, radius, [np.random.rand()*np.pi*2, 0, 0], dock_rotations), mol.FixedSubstrate(substrate_location_5, radius, [np.random.rand()*np.pi*2, 0, 0], dock_rotations), mol.FixedSubstrate(substrate_location_6, radius, [np.random.rand()*np.pi*2, 0, 0], dock_rotations), mol.FixedSubstrate(substrate_location_7, radius, [np.random.rand()*np.pi*2, 0, 0], dock_rotations), mol.FixedSubstrate(substrate_location_8, radius, [np.random.rand()*np.pi*2, 0, 0], dock_rotations), mol.FixedSubstrate(substrate_location_9, radius, [np.random.rand()*np.pi*2, 0, 0], dock_rotations)]

      ligand_orientation = [substrates_start[0].start_rotation[0], np.pi, 0]
      ligand_start = mol.Ligand(ligand_location, radius, ligand_orientation, dock_rotations)
      
      run_reactions(file_name, substrates_start, ligand_start)
