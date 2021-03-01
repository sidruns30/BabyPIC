"""
Particle Pusher - See documentation for implementation details 
"""
from math import sqrt
import numpy as np
from interpolate_fields import interpolate_fields

def push(FieldsObject, FpfeObject, dt=None):
    if not dt:
        dt = FieldsObject.dt
    m = FpfeObject.m
    q = FpfeObject.q
    Ex, Ey, B = interpolate_fields(FpfeObject, FieldsObject)
    E = np.array([Ex,Ey])    
    # Useful variables 
    gamma = FpfeObject.gamma_new
    # 3 velocity 
    u = gamma*FpfeObject.v_new 
    v_minus = u + (q*dt/m/2)*E
    # Invert matrix to find v_plus
    b = np.array([v_minus[0]/dt + q*B/2/m/gamma*v_minus[1], v_minus[1]/dt -
                  q*B/2/m/gamma*v_minus[0]])
    det_A = (1/dt**2 + (q*B/2/m/gamma)**2)
    v_plus = np.array([b[0]/dt + (q*B/2/m/gamma)*b[1], (-q*B/2/m/gamma)*b[0] 
                       + b[1]/dt])/det_A 
    u_new = v_plus + q*dt/2/m*E
    # Convert to regular velocity 
    v_new = u_new/sqrt(1 + u_new[0]**2 + u_new[1]**2)
    FpfeObject.v_old = FpfeObject.v_new
    FpfeObject.v_new = v_new
    FpfeObject.gamma_old = gamma
    FpfeObject.gamma_new = 1/sqrt(1 - v_new[0]**2 - v_new[1]**2)

    