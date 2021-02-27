"""
Define the class 'FPFE' that stores details about the individual 
FPFE cells in the simulation. The dimension of each FPFE cell is 
the same as that of a 'FIELDS' grid cell. The input of the FPFE
are objects GridParam and fpfe_prop
"""
import numpy as np
from numpy.linalg import norm

class FPFE:
    def __init__(self, GridParam, q, m, r_new, v_new):
        self.q = q
        self.m = m
        self.r_old = np.zeros(3)
        self.r_new = r_new
        self.v_old = np.zeros(3)
        self.v_new = v_new
        self.gamma_old = 1/np.sqrt(1 - norm(self.v_old)**2)
        self.gamma_new = 1/np.sqrt(1 - norm(self.v_new)**2)
        self.dx = (GridParam.xmax - GridParam.xmin)/GridParam.nx
        self.dy = (GridParam.ymax - GridParam.ymin)/GridParam.ny
        # Timestep is limited by speed of light taken from fields_class
        self.dt = 0.05/(np.sqrt(1/self.dx**2 + 1/self.dy**2))
    # Update position and velocity using leapfrog integration
    # Update postion before the velocity 
    def update_pos(self):
        self.r_old = self.r_new
        self.r_new = self.r_old + self.dt*self.v_new
        