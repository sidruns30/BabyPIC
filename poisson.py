"""
Solve the Poisson Equation to find the electric fields from
the intial FPFE distribution.
This function is called once at the start of the simulation
"""
import numpy as np

### Deposit charge on the Yee Grid from FPFE positions 
def charge_deposition(FieldsObject, FPFEList):
    print('Deposting Charges...')
    dx = FieldsObject.dx
    dy = FieldsObject.dy
    xmin = FieldsObject.xmin
    ymin = FieldsObject.ymin
    x_edges = FieldsObject.x_edges
    y_edges = FieldsObject.y_edges
    for fpfe in FPFEList:
        # fpfe top right corner cell index
        i = int((fpfe.r_new[0] - xmin + dx/2)//dx)
        j = int((fpfe.r_new[1] - ymin + dy/2)//dy)
        # Calculate area weights 
        Sx2 = (fpfe.r_new[0] + dx/2 - x_edges[i])/dx
        Sx1 = 1 - Sx2
        Sy2 = (fpfe.r_new[1] + dy/2 - y_edges[j])/dy
        Sy1 = 1 - Sy2    
        
        rho = fpfe.q / (dx*dy)
        FieldsObject.rho[i,j] += Sx2*Sy2*rho
        FieldsObject.rho[i-1,j] += Sx1*Sy2*rho
        FieldsObject.rho[i,j-1] += Sx2*Sy1*rho
        FieldsObject.rho[i-1,j-1] += Sx1*Sy1*rho
        
### Calculate Potential from charge densities
def poisson_solver(FieldsObject):
    nx = FieldsObject.nx
    ny = FieldsObject.ny
    dx = FieldsObject.dx
    dy = FieldsObject.dy  
    # Potential defined at cell centers
    V = np.zeros((nx+2,ny+2))
    k = (dx**2*dy**2)/(-2*(dx**2 + dy**2))
    converge = False
    while not converge:
        Vnew = k*(-FieldsObject.rho[1:-1,1:-1] - 
                  ((V[2:,1:-1] + V[:-2,1:-1])/dx**2 + 
                   (V[1:-1,2:] + V[1:-1,:-2])/dy**2))
        # Check if values are within 1%
        converge = np.allclose(Vnew, V[1:-1,1:-1], rtol=1e-3)
        V[1:-1,1:-1] = Vnew
    # Get potenential at cell corners 
    Vcor = (V[:-1,1:] + V[:-1,:-1] + V[1:,:-1] + V[1:, 1:])/4 
    return Vcor

### Calculate the Electric fields from the gradient of the potential
### Also returns the potential 
def grad_v(FieldsObject):
    print('Solving Poisson Equation...')
    dx = FieldsObject.dx
    dy = FieldsObject.dy
    Vcor = poisson_solver(FieldsObject)
    print('Calculating Electric Fields...')
    FieldsObject.Ex_new[1:-1,:] = -(Vcor[1:,:] - Vcor[:-1,:])/dx
    FieldsObject.Ey_new[:,1:-1] = -(Vcor[:,1:] - Vcor[:,:-1])/dy

    
    
        
        
        
        