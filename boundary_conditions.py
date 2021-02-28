"""
Set boundary conditions for the stationary grid and the FPFEs
"""

def set_boundary_condition_fields(FieldsObject, step, bx='periodic', by='periodic'):
    if bx == 'periodic':
        if step == 0:
            FieldsObject.Ex_new[-1, :] = FieldsObject.Ex_new[1, :]
            FieldsObject.Ex_new[0, :] = FieldsObject.Ex_new[-2, :]
            FieldsObject.Ey_new[-1, :] = FieldsObject.Ey_new[1, :]
            FieldsObject.Ey_new[0, :] = FieldsObject.Ey_new[-2, :]
        if step == 1/2:
            FieldsObject.Bz_new[-1, :] = FieldsObject.Bz_new[1, :]            
            FieldsObject.Bz_new[0, :] = FieldsObject.Bz_new[-2, :]
    if by == 'periodic':
        if step == 0:
            FieldsObject.Ex_new[:, -1] = FieldsObject.Ex_new[:, 1]
            FieldsObject.Ex_new[:, 0] = FieldsObject.Ex_new[:, -2]
            FieldsObject.Ey_new[:, -1] = FieldsObject.Ey_new[:, 1]
            FieldsObject.Ey_new[:, 0] = FieldsObject.Ey_new[:, -2]     
        if step == 1/2:
            FieldsObject.Bz_new[:, -1] = FieldsObject.Bz_new[:, 1]            
            FieldsObject.Bz_new[:, 0] = FieldsObject.Bz_new[:, -2]
    if bx == 'conductiing':
        if step == 0:
            FieldsObject.Ex_new[0, :] = 0
            FieldsObject.Ex_new[-1, :] = 0
    if by == 'conducting':
        if step == 0:
            FieldsObject.Ey_new[:, 0] = 0
            FieldsObject.Ey_new[:, -1] = 0

def set_boundary_condition_fpfe(FieldsObject, FPFEObject, step, bx='periodic', by='periodic'):
    xmin = FieldsObject.xmin
    xmax = FieldsObject.xmax
    ymin = FieldsObject.ymin
    ymax = FieldsObject.ymax
    if bx == 'periodic':
        if step == 0: 
            delx = 0
            FPFEObject.r_new[0] = FPFEObject
        if step == 1/2:
            pass
    if by == 'periodic':
        if step == 0:
            pass
        if step == 1/2:
            pass
    if bx == 'reflecting':
        if step == 0: 
            pass
        if step == 1/2:
            pass
    if by == 'reflecting':
        if step == 0:
            pass
        if step == 1/2:
            pass
