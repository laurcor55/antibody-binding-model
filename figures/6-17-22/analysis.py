import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
file_name = 'output_moveback.json'
file = open('output_moveback.json', 'r')
df = pd.json_normalize(json.loads(file.read()))
#df['substrates'].apply(pd.Series)
df_substrates = pd.json_normalize(df['substrates'])
df = pd.concat([df, df_substrates], axis=1)
df = df.drop(columns = ['substrates'])
df.columns = df.columns.map(str)
df_substrate = pd.json_normalize(df['0'])
df = pd.concat([df, df_substrate], axis=1)
df.info()

df['p_bind'] = df['reaction.binding_count']/df['reaction.total_count']
df['location'] = df['location'].apply(np.asarray)
df['ligand.location'] = df['ligand.location'].apply(np.asarray)

df['closest_center_distance'] = np.linalg.norm(df['location'] - df['ligand.location'])

#reactions = pd.read_json(file_name)
#reactions = pd.json_normalize(reactions['reaction'])
#reactions = reactions['reaction'].apply(pd.Series)
n_reactions = len(reactions)
#print(reactions[0]['molecules'][0]['location'])
#print(overall_list[0]['molecules'][0]['n_docks'])
p_bind = []
center_distance = []
n_docks = []
n_substrates = []
for ii in range(n_reactions):
  reaction = reactions[ii]
  p_bind.append(reaction['reaction']['binding_count']/reaction['reaction']['total_count'])
  ligand_location = np.asarray(reaction['substrates'][0]['location'])
  substrate_location = np.asarray(reaction['ligand']['location'])

  center_distance.append( np.linalg.norm(ligand_location - substrate_location))
  
  n_substrates.append(len(reaction['substrates']))
  n_docks.append(reaction['ligand']['n_docks'])

plt.scatter(center_distance, p_bind, c=n_docks)
plt.show()

plt.scatter(center_distance, p_bind, c=n_substrates)
plt.show()