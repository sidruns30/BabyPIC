"""
Define the class 'FIELDS' that holds electric and magnetic fields
along with currents, initial charge densities and grid parameters
Speed of light is assumed to be 1 
"""
import numpy as np

class FIELDS():
    def __init__(self, GridParam):
        '''
        GridParam: Object (Contains cell number, grid limits and 
        boundary conditions for the fields)
        '''
        self.nx = GridParam.nx
        self.ny = GridParam.ny
        self.dx = (GridParam.xmax - GridParam.xmin)/GridParam.nx
        self.dy = (GridParam.ymax - GridParam.ymin)/GridParam.ny
        # FIELDS coordinate limit includes ghost cells 
        self.xmin = GridParam.xmin - self.dx
        self.xmax = GridParam.xmax + self.dx
        self.ymin = GridParam.ymin - self.dy
        self.ymax = GridParam.ymax + self.dy
        # Timestep from Courant condition
        self.dt = 0.05/(np.sqrt(1/self.dx**2 + 1/self.dy**2))
        # Coordinates of cell edges 
        self.x_edges = np.linspace(self.xmin, self.xmax, self.nx+3) 
        self.y_edges = np.linspace(self.ymin, self.ymax, self.ny+3)
        # Coordinates of cell centers 
        self.x_cen = (self.x_edges[1:] + self.x_edges[:-1])/2 
        self.y_cen = (self.y_edges[1:] + self.y_edges[:-1])/2 
        # Set up zero fields, currents 
        self.Bz_old = 1e-10*np.zeros((self.nx + 2, self.ny + 2))
        self.Ex_old = 1e-10*np.zeros((self.nx + 2, self.ny + 1))
        self.Ey_old = 1e-10*np.zeros((self.nx + 1, self.ny + 2))
        self.Bz_new = 1e-10*np.zeros((self.nx + 2, self.ny + 2))
        self.Ex_new = 1e-10*np.zeros((self.nx + 2, self.ny + 1))
        self.Ey_new = 1e-10*np.zeros((self.nx + 1, self.ny + 2))        
        self.Jx = np.zeros((self.nx + 2, self.ny + 1))
        self.Jy = np.zeros((self.nx + 1, self.ny + 2))
        self.rho = np.zeros((self.nx + 2, self.ny + 2))
    # Update fields using leapfrog integration
    # Update Ex and Ey before updating Bz
    def update_ex(self):
        self.Ex_old = self.Ex_new
        # Partial Bz partial y
        t = self.dt/self.dy*(self.Bz_new[1:-1,2:-1] - self.Bz_new[1:-1,1:-2])
        self.Ex_new[1:-1,1:-1] = (self.Ex_old[1:-1,1:-1] + t - 
                                  self.dt*4*np.pi*self.Jx[1:-1,1:-1])
        self.Jx[:,:] = 0
        
    def update_ey(self):
        self.Ey_old = self.Ey_new
        # Partial Bz partial x 
        t = self.dt/self.dx*(self.Bz_new[2:-1,1:-1] - self.Bz_new[1:-2,1:-1])
        self.Ey_new[1:-1,1:-1] = (self.Ey_old[1:-1,1:-1] - t - 
                                  self.dt*4*np.pi*self.Jy[1:-1, 1:-1])
        self.Jy[:,:] = 0
        
    def update_bz(self):
        self.Bz_old = self.Bz_new
        a = 0.125*self.dx/self.dy
        b = 1 - 2*a
        # Partial Ex partial y (t1 + t2 + t3)
        t1 = b*(self.Ex_new[1:-1,1:] - self.Ex_new[1:-1,0:-1])
        t2 = a*(self.Ex_new[2:,1:] - self.Ex_new[2:,0:-1])
        t3 = a*(self.Ex_new[0:-2,1:] - self.Ex_new[0:-2,0:-1])
        # Partial Ey partial x (t4 + t5 + t6)
        t4 = b*(self.Ey_new[1:,1:-1] - self.Ey_new[0:-1,1:-1])
        t5 = a*(self.Ey_new[1:,2:] - self.Ey_new[0:-1,2:])
        t6 = a*(self.Ey_new[1:,0:-2] - self.Ey_new[0:-1,0:-2])
        self.Bz_new[1:-1,1:-1] = self.Bz_old[1:-1,1:-1] + self.dt*((t1 + t2 + t3)/self.dy - 
                                                                   (t4 + t5 + t6)/self.dx)
