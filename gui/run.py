import tkinter as tk
import matplotlib

matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import mpl_toolkits.mplot3d.axes3d as p3

import sys
sys.path.insert(1, '../')
import reaction as reac
import molecules as mol


class BrownianDynamics():
  def __init__(self, root):
    self.root = root
    self.root.geometry("1000x500")
    self.create_default_reaction()
    self.create_input_frame()
    self.create_output_frame()

  def create_default_reaction(self):
    ligand_location = [0, 0, 42]
    ligand_orientation = [0, 0, 0]
    substrate_location = [0, 0, 0]
    substrate_orientation = [0, 0, 0]
    dock_rotations = [[0, 0, 0]]
    radius = 18
    rend = 200
    minimum_binding_docks = 3
    ligand = mol.Ligand(ligand_location, radius, ligand_orientation, dock_rotations)
    substrates = [mol.FixedSubstrate(substrate_location, radius, substrate_orientation, dock_rotations)]
    self.reaction = reac.Reaction(ligand, substrates, rend, minimum_binding_docks)

  def create_input_frame(self):
    self.input_frame = tk.Frame()
    label = tk.Label(self.input_frame, text='Input Parameters')
    label.pack(padx=10, pady=10)
    fig = Figure(figsize=(5, 5), dpi=100)
    
    self.start_view_axis = p3.Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(self.start_view_axis)
    self.reaction.show_molecules_at_time(self.start_view_axis, 0)
    
    self.start_view_canvas = FigureCanvasTkAgg(fig, self.input_frame)
    self.start_view_canvas.draw()
    self.start_view_canvas.get_tk_widget().pack()
    
    self.start_button = tk.Button(self.input_frame, text='Start Reaction', command=self.run_reaction)
    self.start_button.pack()

    self.input_frame.grid(row=0, column=0)

  def run_reaction(self):
    self.reaction.back_to_start()
    self.reaction.progress_reaction()
    self.reaction.show_animation(self.output_figure)
    self.reaction_view_canvas.draw()


  def create_output_frame(self):
    self.output_frame = tk.Frame()
    label = tk.Label(self.output_frame, text='output stuff')
    label.pack(padx=10, pady=10)
    self.output_figure = Figure(figsize=(5, 5), dpi=100)
    
    self.reaction_view_canvas = FigureCanvasTkAgg(self.output_figure, self.output_frame)
    self.reaction_view_canvas.draw()
    self.reaction_view_canvas.get_tk_widget().pack()
    self.output_frame.grid(row=0, column=1)

if __name__ == '__main__':
  root = tk.Tk()

  app = BrownianDynamics(root)
  root.mainloop()