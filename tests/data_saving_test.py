import sys
sys.path.insert(1, '../')

import json
import molecules as mol
import reaction as reac
import numpy as np

def append_output(new_data, file_name):
  with open(file_name, 'r+') as file:
    file_data = json.load(file)
    file_data.append(new_data)
    file.seek(0)
    json.dump(file_data, file, indent=4)

def write_output(output_dict, file_name):
  with open(file_name, "w") as outfile:
    json.dump(output_dict, outfile)

def to_dict(molecules_dict, reaction_dict):
  overall_json = {}
  overall_json['reaction'] = reaction_dict
  overall_json['molecules'] = molecules_dict
  return overall_json

output_list = []
radius = 18
rend = 200
minimum_binding_docks = 3
ligand_location = [0, 0, 39]
substrate_location = [0, 0, 0]
ligand_orientation = [0, np.pi, 0]
substrate_orientation = [0, 0, 0]
n_docks = 8
jj = 0
while jj < 10:
  molecules = [mol.Ligand(ligand_location, radius, ligand_orientation, n_docks), mol.FixedSubstrate(substrate_location, radius, substrate_orientation, n_docks)]
  molecules_dict = [molecule.__dict__() for molecule in molecules]

  reaction = reac.Reaction(molecules, rend, minimum_binding_docks, 'end_at_binding')
  reaction.progress_reaction()
  reaction_dict = reaction.__dict__()

  output_list.append(to_dict(molecules_dict, reaction_dict))
  jj += 1

write_output(output_list, 'overall.json')
#append_output(output_list, 'overall.json')