import numpy as np
from matplotlib import pyplot as plt

class scaling_data:

    def __init__(self, axis_cells=[0,0,0], dim=0, nproc=0, time_sec=0.0, \
        mesh_block=[0,0,0],pos_rseed=False):
        self.axis_cells = np.array(axis_cells)  # Number of cells per axis
        self.dim    = dim                       # Number of dimensions
        self.ncells = np.prod(axis_cells)       # Total number of cells
        self.nproc  = nproc                     # Number of processors used
        self.time   = time_sec                  # Time taken in seconds
        self.pos_rseed = pos_rseed

        if mesh_block==[0,0,0]:
            self.mesh_block = self.axis_cells
        else:
            self.mesh_block = np.array(mesh_block) # Meshblock for parallel processing

    def time_min(self):
        return self.time/60
    
    def time_hr (self):
        return self.time_min()/60

    def show_data(self):
        print("Number of cells per axis : ",self.axis_cells)
        print("Number of dimensions     : ",self.dim)
        print("Number of total cells    : ",self.ncells)
        print("Number of processors     : ",self.nproc)
        print("Time taken in seconds    : ",self.time,"sec")
        print("Meshblock                : ",self.mesh_block)
        print("Is rseed positive?       : ",self.pos_rseed)