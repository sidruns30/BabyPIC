"""
Set initial conditions, Solve Gauss's Law'
Input values: electric fields in the background (0 divergence) t = 0
              positions of the charges t = 0
              magnetic fields in the background (0 divergence) t = 1/2 dt
              velocities of the charges t = 1/2
TIMELINE:
    STEP 0 Solve Poisson equation, calculate and add electric fields from the
           charges to the diveregence free background field
           --------------Main Simulation Start-----------------
    STEP 1: Update particle positions (t=n+1) from velocities(t=n+1/2)
    STEP 2: Deposit currents (t=n+1/2) on to the grid from average of new and
            old positions (t=n,t=n+1) and momenta (t=n+1/2)
    STEP 3: Calculate electric fields (t=n+1) from currents and magnetic 
            field (t=n+1/2)
    STEP 4: Calculate magnetic fields (t=n+3/2) from new electric 
            fields (t=n+1)
    STEP 5: Update the momentum (t=n+3/2) of the particles from time-averaged 
            magnetic fields (t=n+1/2,t=n+3/2) and electric fields (t=n+1)

NOTES:
    - FPFE and particles refer to the same thing 
    - 'push()' is yet to be tested 
    - Equations are written with units of c=1, mu=1 and ep=1
"""

from fields_class import FIELDS
from fpfe_class import FPFE
from poisson import charge_deposition, grad_v
from particle_push import push
from current_deposition import deposit
from interpolate_fields import interpolate_fields
from boundary_conditions import set_boundary_condition_fields
from boundary_conditions import set_boundary_condition_fpfe

import numpy as np
from numpy.linalg import norm 

###---------------INPUT-----------------------
# Set the grid properties 
class GRID:
    def __init__(self):
        self.nx = 200
        self.ny = 200
        self.xmin = -1
        self.xmax = 1
        self.ymin = -1
        self.ymax = 1
GridParam = GRID()

### Set the FPFE properties 
N_FPFE = 1
charges = [1,1]
masses = [1,1]
# Set positions
positions = np.array([[-0,0], [2,0]])
velocities = np.array([[0,0], [0,0]])
FPFE_list = []
for i in range(N_FPFE):
    FPFE_list.append(FPFE(GridParam, q=charges[i], m=masses[i], 
                          r_new=positions[i], v_new=velocities[i]))

###-----------------SIMULATION START---------------- 
# Create the FieldsOject
FieldsObject = FIELDS(GridParam)

'''
STEP 0
'''
# Solve Gausses Law and calculate electric fields
charge_deposition(FieldsObject, FPFE_list)
grad_v(FieldsObject)
# Add divergence free background electric field to existing field
Ex_div_free = -100*np.ones((GridParam.nx+2, GridParam.ny+1))
Ey_div_free = np.zeros((GridParam.nx+1, GridParam.ny+2))

FieldsObject.Ex_new[:,:] += Ex_div_free
FieldsObject.Ey_new[:,:] += Ey_div_free

E1 = []
E2 = []
B1 = []
t = []
### Define one simulation timestep 
def timestep(FieldsObject, FPFE_list):
    step = 0
    # STEPS 1 and 2 - Update particle positions + deposit currents on the grid
    for FpfeObject in FPFE_list:
        FpfeObject.update_pos()
        deposit(FieldsObject, FpfeObject)
    # STEP 3 - Update Electric Fields 
    FieldsObject.update_ex()
    FieldsObject.update_ey()
    set_boundary_condition_fields(FieldsObject, step, 'periodic', 'periodic')
    # STEP 4 - Update Magnetic Fields 
    step = 1/2
    FieldsObject.update_bz()
    set_boundary_condition_fields(FieldsObject, step, 'periodic', 'periodic')
    # STEP 5 - Update particle velocities 
    for FpfeObject in FPFE_list:
        push(FieldsObject, FpfeObject)
        

timesteps = 10000
print('Running Main Simulation...')
FpfeObject = FPFE_list[0]
for i in range(timesteps):
    print('Timestep: ', i+1)
    timestep(FieldsObject, FPFE_list)
    print('Velocity Magnitudes:')
    print([norm(FPFE_list[0].v_new)])#, norm(FPFE_list[1].v_new)])
    t.append(i)
    e1, e2, b1 = interpolate_fields(FPFE_list[0], FieldsObject)
    E1.append(e1)
    E2.append(e2)
    B1.append(b1)
'''
    if (i%100==0):
        
        
        ###----------------PLOTTING------------------------------- 
        
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10,10))
        X,Y = np.meshgrid(FieldsObject.x_edges, FieldsObject.y_edges)
        for FpfeObject in FPFE_list:
            plt.plot(FpfeObject.r_new[0], FpfeObject.r_new[1], 'go', alpha=1)
        #plt.pcolormesh(X,Y,FieldsObject.rho.T)
        plt.pcolormesh(X[1:-1,1:-1],Y[1:-1,1:-1],FieldsObject.Bz_old[1:-1,1:-1].T, 
                       cmap='seismic', vmin=-5, vmax=5)
        plt.colorbar()
        
        Ex = (FieldsObject.Ex_old[1:,:] + FieldsObject.Ex_old[:-1,:])/2
        Ey = (FieldsObject.Ey_old[:,1:] + FieldsObject.Ey_old[:,:-1])/2
        plt.streamplot(X[1:-1,1:-1],Y[1:-1,1:-1],Ex.T,Ey.T, color='k')
        plt.xlim((GridParam.xmin, GridParam.xmax))
        plt.ylim((GridParam.ymin, GridParam.ymax))
        plt.title('$B_z$', fontsize=20)
        plt.savefig('plots/fig'+str(i//100)+'.png')
        plt.close()

E1 = np.array(E1)
E2 = np.array(E2)
'''
import matplotlib.pyplot as plt

plt.plot(t,E1,t,E2,t,)

plt.plot(t,E1,label='Ex')
plt.plot(t,E2,label='Ey')
plt.plot(t,B1,label='Bz')
plt.legend()





