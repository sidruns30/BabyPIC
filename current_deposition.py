'''
CURRENT INTERPOLATION FROM FPFEs ON TO THE GRID
'''
### Break displacement if FPFE crosses a cell (Umeda, et.al)
def deposit(FieldsObject, FpfeObject):
    dx = FieldsObject.dx
    dy = FieldsObject.dy
    xmin = FieldsObject.xmin
    ymin = FieldsObject.ymin
    x_old = FpfeObject.r_old[0]
    y_old = FpfeObject.r_old[1]
    x_new = FpfeObject.r_new[0]
    y_new = FpfeObject.r_new[1]
    q = FpfeObject.q
    # Cell indices for initial and final positions of fpfe
    i1 = int((x_old - xmin)//dx)
    i2 = int((x_new - xmin)//dx)
    j1 = int((y_old - ymin)//dy)
    j2 = int((y_new - ymin)//dy)
    # Intermediate positions 
    xr = (x_old + x_new)/2 if i1 == i2 else max(i1*dx, i2*dx) + xmin
    yr = (y_old + y_new)/2 if j1 == j2 else max(j1*dy, j2*dy) + ymin
    dA = dx * dy
    dt = FieldsObject.dt
    # Calculate weights
    Wx1 = (x_old + xr)/2/dx - xmin/dx - i1
    Wx2 = (xr + x_new)/2/dx - xmin/dx - i2
    Wy1 = (y_old + yr)/2/dy - ymin/dy - j1
    Wy2 = (yr + y_new)/2/dy - ymin/dy - j2
    # Calculate the charge fluxes
    Fx1 = q*(xr - x_old)/dt
    Fx2 = q*(x_new - xr)/dt
    Fy1 = q*(yr - y_old)/dt
    Fy2 = q*(y_new - yr)/dt
    # Deposit currents 
    FieldsObject.Jx[i1,j1-1] += Fx1*(1 - Wy1)/dA
    FieldsObject.Jx[i1,j1] += Fx1*Wy1/dA
    FieldsObject.Jy[i1-1,j1] += Fy1*(1 - Wx1)/dA
    FieldsObject.Jy[i1,j1] += Fy1*Wx1/dA
    
    FieldsObject.Jx[i2,j2-1] += Fx2*(1 - Wy2)/dA
    FieldsObject.Jx[i2,j2] += Fx2*Wy2/dA
    FieldsObject.Jy[i2-1,j2] += Fy2*(1 - Wx2)/dA
    FieldsObject.Jy[i2,j2] += Fy2*Wx2/dA



