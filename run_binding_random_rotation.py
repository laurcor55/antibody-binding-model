#kernprof -l run.py
#lineprofilergui -l run.py.lprof
import reaction as reac
import molecules as mol
import numpy as np
import time
import json
import copy

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

def run_reactions_random_rotation(file_name, substrates, ligand):
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
  
  while binding_count < 20:
  #while kk < total_count:
    reaction.reset_molecules_to_rotation()
    reaction.progress_reaction()
    reaction.show_animation()
    if reaction.reaction_status == 'binding':
      binding_count += 1
      if (reaction.substrates[0].locked == True):
        center_binding_count += 1
    if kk % 10 ==0:
      k_on = calculate_kon(binding_count, kk+1, ligand, substrates, rstart, rend)
      print(str(kk) + ' of ' + str(total_count) + ', ' + str(binding_count) + ' bound (' + reac.scientific(binding_count/(kk+1)) + '), ' +  k_on  + ' /M/s', end="\r")
    kk += 1
  total_count = kk
  center_binding_proportion = center_binding_count/total_count

  reaction_dict = {'total_count': total_count, 'binding_count': binding_count, 'minimum_binding_docks': minimum_binding_docks, 'center_binding_proportion': center_binding_proportion}
  output_list = to_dict(ligand_dict, substrates_dict, reaction_dict)
  
  print(end='\n')
  print(binding_count)
  k_on = calculate_kon(binding_count, total_count, ligand, substrates, rstart, rend)
  print(k_on)
  append_output(output_list, file_name)
  #write_output(output_list, file_name)

  
file_name = 'figures/6-23-22/output.json'

radius = 18 # angstroms
z_spacing = 40
substrate_spacing = 36

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

ligand_orientation = [0, np.pi, 0]
substrate_orientation = [0, 0, 0]
n_docks = 4

substrates_start = [mol.FixedSubstrate(substrate_location, radius, substrate_orientation, n_docks)]#, mol.FixedSubstrate(substrate_location_2, radius, substrate_orientation, n_docks), mol.FixedSubstrate(substrate_location_3, radius, substrate_orientation, n_docks), mol.FixedSubstrate(substrate_location_4, radius, substrate_orientation, n_docks), mol.FixedSubstrate(substrate_location_5, radius, substrate_orientation, n_docks), mol.FixedSubstrate(substrate_location_6, radius, substrate_orientation, n_docks), mol.FixedSubstrate(substrate_location_7, radius, substrate_orientation, n_docks), mol.FixedSubstrate(substrate_location_8, radius, substrate_orientation, n_docks), mol.FixedSubstrate(substrate_location_9, radius, substrate_orientation, n_docks)]
ligand_start = mol.Ligand(ligand_location, radius, ligand_orientation, n_docks)
    
run_reactions_random_rotation(file_name, substrates_start, ligand_start)
