import tkinter as tk
from tkinter import ttk
import matplotlib

matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import mpl_toolkits.mplot3d.axes3d as p3

import sys
import reaction as reac
import molecules as mol
import numpy as np
import time
import copy

class BrownianDynamics():
  def __init__(self, root):
    self.root = root
    self.root.geometry("1400x900")
    self.create_default_reaction()
    self.layout_app()

  def layout_app(self):
    self.root.columnconfigure(0, weight=1)
    self.root.columnconfigure(1, weight=1)

    style = ttk.Style()
    #style.configure('TFrame', background='green')
    #style.configure('TFrame', background='gray')

    self.create_input_frame()
    self.create_progress_frame()
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
    self.total_reactions_run = 0
    self.binding_count = 0

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
    self.input_frame.grid(row=0, column=0, rowspan=6)

    self.preview_button = tk.Button(self.input_frame, text='Preview Settings', command=self.apply_settings)
    self.preview_button.pack()

  def create_progress_frame(self):
    self.progress_frame = ttk.Frame(padding=(20, 20))
    tk.Label(self.progress_frame, text='Reaction Progress').pack(padx=10, pady=10)
    
    self.start_button = tk.Button(self.progress_frame, text='Start Reaction', command=self.run_reactions)
    self.start_button.pack()

    self.stop_button = tk.Button(self.progress_frame, text='Stop Reaction', command=self.stop_reactions)
    self.stop_button.pack()

    self.progress_bar = ttk.Progressbar(self.progress_frame, orient='horizontal', mode='determinate', length=280)
    self.progress_bar.pack()
    self.progress_frame.grid(row=6, column=0, sticky='n', rowspan=1)
      
  def create_output_frame(self):
    self.output_frame = ttk.Frame(padding=(20, 20))
    tk.Label(self.output_frame, text='Output').pack(padx=10, pady=10)
    
    self.create_output_figure()
    self.output_frame.grid(row=0, column=1, sticky='n', rowspan=7)

    self.create_output_text()
  
  
  def create_input_values(self):
    self.input_tab_control = ttk.Notebook(self.input_frame)
    self.input_tab_control.pack(pady=10, expand=True)
    
    self.create_ligand_input()
    self.create_substrates_input()
    self.create_reaction_count_input()
  
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

  def create_reaction_count_input(self):
    reaction_count_tab_frame = ttk.Frame(self.input_tab_control)#, width=400, height=280)
    reaction_count_tab_frame.pack()#fill='both', expand=True)
    self.input_tab_control.add(reaction_count_tab_frame, text='Reaction Count')
    
    tk.Label(reaction_count_tab_frame, text="Total Reactions").grid(row=0, column=0, columnspan=2)
    self.reaction_count_sv = tk.StringVar(self.root, value='1000')
    reaction_count_entry = tk.Entry(reaction_count_tab_frame, textvariable=self.reaction_count_sv, width=10)
    reaction_count_entry.grid(row=0, column=3)
   
  def apply_settings(self):
    ligand_location = [0, 0, 42]
    ligand_orientation = [0, 0, 0]
    for ii in range(3):
      ligand_location[ii] = float(self.ligand_location_sv[ii].get())
      ligand_orientation[ii] = float(self.ligand_rotation_sv[ii].get())*0.0174533
    self.reaction_count = int(self.reaction_count_sv.get())
    dock_rotations = [[0, 0, 0]]
    radius = 18
    rend = 200
    minimum_binding_docks = 3
    self.ligand = mol.Ligand(ligand_location, radius, ligand_orientation, dock_rotations)

    substrate_spacing = float(self.substrate_spacing_sv.get())
    n_substrates = int(self.n_substrates_sv.get())
    self.substrates = self.create_substrates(n_substrates, substrate_spacing)

    self.reaction = reac.Reaction(self.ligand, self.substrates, rend, minimum_binding_docks)
    self.reaction.back_to_start()
    self.reaction.show_molecules_at_time(self.start_view_axis, 0)
  
  def stop_reactions(self):
    global running
    running = False

  def end_reaction_updates(self):
    self.update_output_text()
    if not hasattr(self, 'reaction_copy'):
      self.set_animation()
    
  def run_reactions(self):
    self.apply_settings()
    global running
    running = True
    if hasattr(self, 'reaction_copy'):
      del self.reaction_copy
    total_count = self.reaction_count
    binding_count = 0
    seed = int(time.time())
    np.random.seed(seed)
    kk = 0
    center_binding_count = 0
    #while binding_count < 40:
    while kk < total_count and running:
      self.reaction.back_to_start()
      self.reaction.progress_reaction()
      if self.reaction.reaction_status == 'binding':
        binding_count += 1
        if binding_count == 1:
          self.set_animation()
          self.reaction.no_save_preview()
        if (self.reaction.substrates[0].locked == True):
          center_binding_count += 1
      if kk % 10 ==0:
        print(str(kk) + ' of ' + str(total_count) + ', ' + str(binding_count) + ' bound (' + reac.scientific(binding_count/(kk+1)) + '), ', end="\r")
        self.progress_bar['value'] = kk/total_count*100
        self.root.update()
      kk += 1
      self.total_reactions_run = kk
      self.binding_count = binding_count
    print(str(kk) + ' of ' + str(total_count) + ', ' + str(binding_count) + ' bound (' + reac.scientific(binding_count/(kk)) + '), ', end="\r")
    self.end_reaction_updates()

  def set_animation(self):
    self.reaction_copy = copy.deepcopy(self.reaction)
    self.output_figure.clf()
    self.reaction_copy.show_animation(self.output_figure)
    self.reaction_view_canvas.draw()

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

  def create_output_figure(self):
    self.output_figure = Figure(figsize=(5, 6), dpi=100)
    self.reaction_view_axis = p3.Axes3D(self.output_figure, auto_add_to_figure=False)
    self.output_figure.add_axes(self.reaction_view_axis)
    
    output_figure_frame = tk.Frame(self.output_frame)
    self.reaction_view_canvas = FigureCanvasTkAgg(self.output_figure, output_figure_frame)
    self.reaction_view_canvas.draw()
    self.reaction_view_canvas.get_tk_widget().pack()
    output_figure_frame.pack()
  
  def create_output_text(self):
    output_text_frame = tk.Frame(self.output_frame)
    
    output_text_frame.pack()
    self.output_text = tk.StringVar()
    text = "Total reactions run: " + str(self.total_reactions_run) + '\nTotal binding: ' + str(self.binding_count) + '\nBinding proportion: '
    self.output_text.set(text)
    tk.Label(output_text_frame, textvariable=self.output_text).grid(row=0, column=0, columnspan=2)
  
  def update_output_text(self):
    text = "Total reactions run: " + str(self.total_reactions_run) + '\nTotal binding: ' + str(self.binding_count) + '\nBinding proportion: ' + reac.scientific(self.binding_count/self.total_reactions_run) 
    self.output_text.set(text)

if __name__ == '__main__':
  root = tk.Tk()

  app = BrownianDynamics(root)
  root.mainloop()