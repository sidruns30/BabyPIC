"""
Set initial conditions, Solve Gauss's Law'
"""
from fields_class import FIELDS
from fpfe_class import FPFE
from poisson import charge_deposition, grad_v
from particle_push import push
from current_deposition import deposit

import numpy as np

# Set the grid properties 
class GRID:
    def __init__(self):
        self.nx = 200
        self.ny = 200
        self.xmin = -2
        self.xmax = 2
        self.ymin = -2
        self.ymax = 2
GridParam = GRID()

# Set the FPFE properties 
N_FPFE = 2
charges = [-0.5, 0.5]
masses = np.ones(N_FPFE)
positions = [np.array([-1,0]), np.array([1,0])]
velocities = [np.array([0,0]), np.array([0,0])]
# Initial values for the FPFEs
FPFE_list = []
for i in range(N_FPFE):
    FPFE_list.append(FPFE(GridParam, q=charges[i], m=masses[i], 
                          r_old=positions[i], v_old=velocities[i]))

# Create the FieldsOject
FieldsObject = FIELDS(GridParam)
# Solve Gausses Law
charge_deposition(FieldsObject, FPFE_list)
V = grad_v(FieldsObject)
# Update the Magnetic Field
FieldsObject.update_bz()
# Deposit Current 
print('Depositing Current...')
for FpfeObject in FPFE_list:
    deposit(FieldsObject, FpfeObject)
# Update the Electric field
print('Updating Electric Fields...')
FieldsObject.update_ex()
FieldsObject.update_ey()
# Push the particle position and momentum
print('Updating the particle momentum...') 
for FpfeObject in FPFE_list:
    push(FieldsObject, FpfeObject)


### Calculate the remaining timesteps 
timesteps = 0
for i in range(timesteps):
    print('Timestep = ', i)
    print(FpfeObject.r_new, FpfeObject.r_old)
    # Update particle position and deposit current 
    for FpfeObject in FPFE_list:
        FpfeObject.update_pos()
        deposit(FieldsObject, FpfeObject)
    # Calculate electric fields 
    FieldsObject.update_ex()
    FieldsObject.update_ey()
    # Remove the currents 
    FieldsObject.Jx[:,:] = 0
    FieldsObject.Jy[:,:] = 0
    # Calculate Magnetic fields 
    FieldsObject.update_bz()
    # Update Momentum
    for FpfeObject in FPFE_list:
        push(FieldsObject, FpfeObject)

    


### Plotting 
import matplotlib.pyplot as plt
X,Y = np.meshgrid(FieldsObject.x_edges, FieldsObject.y_edges)
for FpfeObject in FPFE_list:
    print(FpfeObject.v_new)
    plt.plot(FpfeObject.r_old[0], FpfeObject.r_old[1], 'o')
#plt.pcolormesh(X,Y,FieldsObject.rho.T)
plt.pcolormesh(X[1:-1,1:-1],Y[1:-1,1:-1],FieldsObject.Bz_old[1:-1,1:-1].T, cmap='seismic')
plt.colorbar()

Ex = (FieldsObject.Ex_old[1:,:] + FieldsObject.Ex_old[:-1,:])/2
Ey = (FieldsObject.Ey_old[:,1:] + FieldsObject.Ey_old[:,:-1])/2
plt.streamplot(X[1:-1,1:-1],Y[1:-1,1:-1],Ex.T,Ey.T, color='k')
plt.show()
