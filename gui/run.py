import tkinter as tk
from tkinter import ttk
import matplotlib
from pyparsing import col

matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import mpl_toolkits.mplot3d.axes3d as p3

import sys
sys.path.insert(1, '../')
import reaction as reac
import molecules as mol
import numpy as np

class BrownianDynamics():
  def __init__(self, root):
    self.root = root
    self.root.geometry("1000x500")
    self.create_default_reaction()
    self.layout_app()

  def layout_app(self):
    self.root.columnconfigure(0, weight=1)
    self.root.columnconfigure(1, weight=1)

    style = ttk.Style()
    style.configure('TFrame', background='green')
    self.create_input_frame()
    self.create_output_frame()

  def create_default_reaction(self):
    ligand_location = [0, 0, 42]
    ligand_orientation = [0, 0, 0]
    
    dock_rotations = [[0, 0, 0]]
    radius = 18
    rend = 200
    minimum_binding_docks = 3
    n_substrates = 9
    substrate_spacing = 36
    ligand = mol.Ligand(ligand_location, radius, ligand_orientation, dock_rotations)
    substrates = self.create_substrates(n_substrates, substrate_spacing)
    self.reaction = reac.Reaction(ligand, substrates, rend, minimum_binding_docks)

  def create_substrates(self, n_substrates, substrate_spacing):
    substrate_location = [0, 0, 0]
    substrate_orientation = [0, 0, 0]
    radius = 18
    dock_rotations = [[0, 0, 0]]#, [0, np.pi/4, 0]]
    substrate_location_2 = [0, substrate_spacing, 0]
    substrate_location_3 = [0, -substrate_spacing, 0]
    substrate_location_4 = [substrate_spacing, 0, 0]
    substrate_location_5 = [-substrate_spacing, 0, 0]
    substrate_location_6 = [substrate_spacing, substrate_spacing, 0]
    substrate_location_7 = [-substrate_spacing, -substrate_spacing, 0]
    substrate_location_8 = [-substrate_spacing, substrate_spacing, 0]
    substrate_location_9 = [substrate_spacing, -substrate_spacing, 0]
    substrates = [mol.FixedSubstrate(substrate_location, radius, substrate_orientation, dock_rotations), mol.FixedSubstrate(substrate_location_2, radius, substrate_orientation, dock_rotations), mol.FixedSubstrate(substrate_location_3, radius, substrate_orientation, dock_rotations), mol.FixedSubstrate(substrate_location_4, radius, substrate_orientation, dock_rotations), mol.FixedSubstrate(substrate_location_5, radius, substrate_orientation, dock_rotations), mol.FixedSubstrate(substrate_location_6, radius, substrate_orientation, dock_rotations), mol.FixedSubstrate(substrate_location_7, radius, substrate_orientation, dock_rotations), mol.FixedSubstrate(substrate_location_8, radius, substrate_orientation, dock_rotations), mol.FixedSubstrate(substrate_location_9, radius, substrate_orientation, dock_rotations)]
    substrates = substrates[:n_substrates]
    return substrates

  def create_input_frame(self):
    self.input_frame = ttk.Frame(padding=(20, 20))
    tk.Label(self.input_frame, text='Input Parameters').pack(padx=10, pady=10)
    
    self.create_input_figure()
    self.create_input_values()
    self.input_frame.grid(row=0, column=0)
    
    self.start_button = tk.Button(self.input_frame, text='Start Reaction', command=self.run_reaction)
    self.start_button.pack()

    self.preview_button = tk.Button(self.input_frame, text='Preview Settings', command=self.apply_settings)
    self.preview_button.pack()
      
  def create_output_frame(self):
    self.output_frame = ttk.Frame(padding=(20, 20))
    tk.Label(self.output_frame, text='Output').pack(padx=10, pady=10)
    
    self.create_output_figure()
    self.output_frame.grid(row=0, column=1, sticky='n')

  
  def create_input_values(self):
    self.input_tab_control = ttk.Notebook(self.input_frame)
    self.input_tab_control.pack(pady=10, expand=True)
    
    self.create_ligand_input()
    self.create_substrates_input()

  
  def create_ligand_input(self):
    ligand_tab_frame = ttk.Frame(self.input_tab_control)#, width=400, height=280)
    ligand_tab_frame.pack()#fill='both', expand=True)
    self.input_tab_control.add(ligand_tab_frame, text='Ligand Input')
    
    tk.Label(ligand_tab_frame, text="Location (x, y, z)").grid(row=0, column=0, columnspan=2)
    self.ligand_location_sv = []
    ligand_location_entry =[]
    for ii in range(3):
      self.ligand_location_sv.append(tk.StringVar(self.root, value=str(self.reaction.ligand.start_location[ii])))
      ligand_location_entry.append(tk.Entry(ligand_tab_frame, textvariable=self.ligand_location_sv[ii], width=10))
      ligand_location_entry[ii].grid(row=0, column=ii+3)
    tk.Label(ligand_tab_frame, text="Rotation (deg, deg, deg)").grid(row=1, column=0, columnspan=2)
    self.ligand_rotation_sv = []
    ligand_rotation_entry =[]
    for ii in range(3):
      self.ligand_rotation_sv.append(tk.StringVar(self.root, value=str(self.reaction.ligand.start_rotation[ii]/0.0174533)))
      ligand_rotation_entry.append(tk.Entry(ligand_tab_frame, textvariable=self.ligand_rotation_sv[ii], width=10))
      ligand_rotation_entry[ii].grid(row=1, column=ii+3)

  def create_substrates_input(self):
    substrates_tab_frame = ttk.Frame(self.input_tab_control)#, width=400, height=280)
    substrates_tab_frame.pack()#fill='both', expand=True)
    self.input_tab_control.add(substrates_tab_frame, text='Substrates Input')
    
    tk.Label(substrates_tab_frame, text="Number of Substrates").grid(row=0, column=0)
    self.n_substrates_sv = tk.StringVar(self.root, value=str(len(self.reaction.substrates)))
    n_substrates_entry = tk.Entry(substrates_tab_frame, textvariable=self.n_substrates_sv)
    n_substrates_entry.grid(row=0, column=1)

    tk.Label(substrates_tab_frame, text="Substrate Spacing").grid(row=1, column=0)
    self.substrate_spacing_sv = tk.StringVar(self.root, value=str(np.sum(np.abs(self.reaction.substrates[0].location - self.reaction.substrates[1].location))))
    substrate_spacing_entry = tk.Entry(substrates_tab_frame, textvariable=self.substrate_spacing_sv)
    substrate_spacing_entry.grid(row=1, column=1)
   
  def apply_settings(self):
    ligand_location = [0, 0, 42]
    ligand_orientation = [0, 0, 0]
    for ii in range(3):
      ligand_location[ii] = float(self.ligand_location_sv[ii].get())
      ligand_orientation[ii] = float(self.ligand_rotation_sv[ii].get())*0.0174533
    
    dock_rotations = [[0, 0, 0]]
    radius = 18
    rend = 200
    minimum_binding_docks = 3
    ligand = mol.Ligand(ligand_location, radius, ligand_orientation, dock_rotations)

    substrate_spacing = float(self.substrate_spacing_sv.get())
    n_substrates = int(self.n_substrates_sv.get())
    substrates = self.create_substrates(n_substrates, substrate_spacing)

    self.reaction = reac.Reaction(ligand, substrates, rend, minimum_binding_docks)
    self.reaction.back_to_start()
    self.reaction.show_molecules_at_time(self.start_view_axis, 0)

  def create_input_figure(self):
    fig = Figure(figsize=(5, 5), dpi=100)
    self.start_view_axis = p3.Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(self.start_view_axis)
    self.reaction.show_molecules_at_time(self.start_view_axis, 0)
    
    input_figure_frame = tk.Frame(self.input_frame)
    self.start_view_canvas = FigureCanvasTkAgg(fig, input_figure_frame)
    self.start_view_canvas.draw()
    self.start_view_canvas.get_tk_widget().pack()
    input_figure_frame.pack()

  def run_reaction(self):
    self.reaction.back_to_start()
    self.reaction.progress_reaction()
    self.reaction.show_animation(self.output_figure)
    self.reaction_view_canvas.draw()

  def create_output_figure(self):
    self.output_figure = Figure(figsize=(5, 6), dpi=100)
    self.reaction_view_axis = p3.Axes3D(self.output_figure, auto_add_to_figure=False)
    self.output_figure.add_axes(self.reaction_view_axis)
    
    output_figure_frame = tk.Frame(self.output_frame)
    self.reaction_view_canvas = FigureCanvasTkAgg(self.output_figure, output_figure_frame)
    self.reaction_view_canvas.draw()
    self.reaction_view_canvas.get_tk_widget().pack()
    output_figure_frame.pack()

if __name__ == '__main__':
  root = tk.Tk()

  app = BrownianDynamics(root)
  root.mainloop()