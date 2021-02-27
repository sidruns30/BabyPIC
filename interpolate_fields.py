"""
Interpolate the value of fields at 'r_old' at time n+1
"""
def interpolate_fields(FpfeObject, FieldsObject):
    x_fpfe, y_fpfe = FpfeObject.r_old[0], FpfeObject.r_old[1]
    dx = FieldsObject.dx
    dy = FieldsObject.dy
    xmin = FieldsObject.xmin
    ymin = FieldsObject.ymin
    x_edges = FieldsObject.x_edges
    y_edges = FieldsObject.y_edges
    x_cen = FieldsObject.x_cen
    y_cen = FieldsObject.y_cen
    # Calculate indices and weights for cell centered values (Bz, Ex[:,], Ey[,:])
    ic = int((x_fpfe + dx/2 - xmin)//dx)
    jc = int((y_fpfe + dy/2 - ymin)//dy)
    Ci2 = (x_fpfe + dx/2 - x_edges[ic])/dx
    Ci1 = 1 - Ci2
    Cj2 = (y_fpfe + dy/2 - y_edges[jc])/dy
    Cj1 = 1 - Cj2
    # Calculate indices and weights for cell edge values (Ex[,:]. Ey[:,])
    ie = int((x_fpfe - xmin)//dx)
    je = int((y_fpfe - ymin)//dy)
    Ei2 = (x_fpfe + dx/2 - x_cen[ie])/dx
    Ei1 = 1 - Ei2
    Ej2 = (y_fpfe + dy/2 - y_cen[je])/dy
    Ej1 = 1 - Ej2
    # Average the magnetic field at time t = n+1
    B_avg = (FieldsObject.Bz_new + FieldsObject.Bz_old)/2
    Ex = FieldsObject.Ex_new
    Ey = FieldsObject.Ey_new
    # Estimate the magnetic field at the fpfe location
    B_fpfe = (Ci2*Cj2*B_avg[ic,jc] + Ci1*Cj2*B_avg[ic-1,jc] + 
              Ci2*Cj1*B_avg[ic,jc-1] + Ci1*Cj1*B_avg[ic-1,jc-1])
    # Estimate Ex at the fpfe location
    Ex_fpfe = (Ci2*Ej2*Ex[ic,je] + Ci1*Ej2*Ex[ic-1,je] + 
               Ci2*Ej1*Ex[ic,je-1] + Ci1*Ej1*Ex[ic-1,je-1]) 
    # Estimate Ey at the fpfe location
    Ey_fpfe = (Ei2*Cj2*Ey[ie,jc] + Ei1*Cj2*Ey[ie-1,jc] + 
               Ei2*Cj1*Ey[ie,jc-1] + Ei1*Cj1*Ey[ie-1,jc-1])
    return Ex_fpfe, Ey_fpfe, B_fpfe
    
