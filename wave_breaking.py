# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 17:46:00 2016

@author: jdorvinen
"""

# Spectral wave propogation from offshore to nearshore



# Propogation of waves via linear wave theory from the nearshore to the wave
# breaking point. Assume parallel bathymetric contours and cross-shore profile
# of Kriebel and Dean.

# Define cross-shore profile geometry
def cross_shore_distance(h, h_T, m, A):
    '''Returns offshore distance x corresponding to given depth h and other
    parameters controlling the beach profile.
    
    Inputs:
        h   = Water depth at x (meters)
        h_T = Depth at which the linear slope is tangent to the concave profile
        m   = Linear beach-face slope
        A   = 2.25*((w**2)/g)**(1/3), from Eq. 15 of Kriebel & Dean (1993).
              A parameter governing profile steepness, valid for sand where 
              0.1mm < d_50 < 0.4mm
    
    Outputs:
        x = offshore distance (m), corresponding to depth h for defined profile
    '''
    x_0 = h_T/(3*m)
        
    if h < h_T:
        x = h/m
    elif h >= h_T:
        x = x_0 + (h/A)**(3/2)
    return x

def linear_propogation(H, T, theta, x_init, m, A, gamma=0.78):
    '''Returns breaking wave height and depth based on an initial nearshore
    wave condition defined by a wave height, period, and angle originating at a
    given distance from the shoreline. Based on linear wave theory and the 
    assumption of shore-parallel dept contours.
    
    Inputs:
        H = Nearshore wave height (m)
        T = Nearshore wave period (T)
        theta = Nearshore wave angle from shore perpendicular (deg) 
                -80 < theta < 80
        x_init = Distance from shoreline where wave originate (m)
        m   = Linear beach-face slope
        A   = 2.25*((w**2)/g)**(1/3), from Eq. 15 of Kriebel & Dean (1993).
              A parameter governing profile steepness, valid for sand where 
              0.1mm < d_50 < 0.4mm
        gamma = Breaker index, usually taken to be 0.78-1.0.
             
    Outputs:
        H_b = Breaking wave height (m)
        h_b = Breaking wave water depth (m)
        x_b = Breaker distance offshore (m)
    '''
    m = m
    A = A 
    gamma = gamma
    
    h_T = (4/9)*(A**3/m**2) # Depth at which the linear slope is tangent to the concave profile
    x_0 = h_T/(3*m)





# Wave set-up calculated from breaking wave height and gamma
def wave_setup(H_b, gamma=0.78):
    '''Returns wave set-up corresponding to a given breaking wave height and 
    breaker index. Eq. 29 in Callaghan et al. (2008), taken originally from 
    Dean and Dalrymple (1991). 
    
    Inputs:
        H_b = Breaking wave height (m)
        gamma = Breaker index, usually taken to be 0.78-1.0.
    
    Outputs:
        nu_max = Maximum wave set-up
    '''
    nu_max = (40-3*gamma**2)*(gamma*H_b)/128
    return nu_max
