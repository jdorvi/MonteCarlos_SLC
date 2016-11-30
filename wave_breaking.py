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
def cross_shore_distance(h, m, A):
    '''Returns offshore distance x corresponding to given depth h and other
    parameters controlling the beach profile. Based on the equilibrium beach
    profile shown by Dean (1976, 1977).
    
    Inputs:
        h   = Water depth at x (meters)
        m   = Linear beach-face slope
        A   = 2.25*((w**2)/g)**(1/3), from Eq. 15 of Kriebel & Dean (1993).
              A parameter governing profile steepness, valid for sand where 
              0.1mm < d_50 < 0.4mm
    
    Outputs:
        x = offshore distance (m), corresponding to depth h for defined profile
    '''
    h_T = (4/9)*(A**3/m**2) # Depth at which the linear slope is tangent to the 
                            # concave profile
    x_0 = h_T/(3*m)
        
    if h < h_T:
        x = h/m
    elif h >= h_T:
        x = x_0 + (h/A)**(3/2)
    return x

def input2deep(H_i, T_i, theta_i, h_i):
    '''Calculates equivalent deep water, offshore, wave conditions given 
    initial wave conditions at a specified depth. Based on Snell's Law and the 
    assumption of shore-parallel depth contours. 
    
    Inputs:
        REQUIRED
        H_i = Initial wave height (m)
        h_i = Initial water depth (m)
        T_i = Initial wave period (sec)
        theta_i = Initial wave angle (deg)
        h_i = Initial water depth (m)
        
    Outputs:
        H_o = Nearshore wave height (m)
        T_o = Nearshore wave period (sec)
        theta_o = Neashore wave angle (deg)
        L_o = Offshore wave length (m)
        C_o = Offshoer wave celerity (m/s)
    '''
    # Required Constants
    from numpy import pi, tanh, arcsin, sin, cos, isnan
    g = 9.8066 #(m/s/s) - Gravitational acceleration
    
    # Convert input angle to radians
    theta_i = theta_i*(2*pi)/360
    
    # Calculate initial wave length, assuming deep water conditions 
    # (ie d/L>0.5)    
    L_i_deep = g*T_i**2/(2*pi)       # Wave length (m)
    
    # Assuming shallow water conditions
    L_i_shallow = (T_i*(g*h_i)**0.5) # Wave length (m)
    
    if h_i/L_i_deep > 0.5:
        # Deep water conditions already at buoy, no need to calculate 
        # equivalent offshore wave parameters
        H_o = H_i
        T_o = T_i
        theta_o = theta_i * 360/(2*pi)
        L_o = L_i_deep
        C_o = ((g*L_o)/(2*pi))**0.5
        return [H_o, T_o, theta_o, L_o, C_o]

    elif h_i/L_i_shallow < 0.05:
        # Shallow water conditions
        C_i = (g*h_i)**0.5
        L_i = L_i_shallow
    
    else:
        # Intermediate depth wave conditions...this is more complicated
        from scipy.optimize import brentq
        
        # Celerity and wave length for intermediate wave conditions
        # Source: Robert M. Sorensen, "Basic Wave Mechanics for Coastal and 
        # Ocean Engineers" (John Wiley & Sons, 1993), Chapter 2.
        def find_L_i(L, T, h):
            # Equation for the length of an intermediate depth wave needs to be
            # solved iteratively
            zero = L - (g*T**2)/(2*pi) * tanh(2*pi*h/L)
            return zero
        L_i = brentq(find_L_i, 10, max(L_i_deep, L_i_shallow), args=(T_i, h_i))
        C_i = ((g*L_i)/(2*pi) * tanh(2*pi*h_i/L_i))**0.5
    
    # Converting to deepwater equivalent wave condition, therefore wave length
    # and wave group celerity are known simply from the wave period.
    L_o = L_i_deep
    C_o = ((g*L_o)/(2*pi))**0.5

    # Calculate shoaling and refraction coefficients following the methodology
    # presented in 'Water Wave Mechanics for Engineers and Scientists' by
    # Dean & Dalrymple (c) 2000
    K_s = (C_i/C_o)**0.5 # Shoaling coefficient. Eq. 4.116
    
    # Find the offshore equivalent wave angle
    theta_o = arcsin(sin(theta_i)*(C_o/C_i)) # Snell's law, Eq. 4.109
    K_r = (cos(theta_i)/cos(theta_o))**0.5 # Refraction coefficient. Eq. 4.118
    
    # Apply shoaling and refraction coefficients to find new wave height
    H_o = H_i*K_s*K_r # Eq. 4.117
    
    # Wave period is assumed constant
    T_o = T_i
    
    # Check for oversteep waves
    if H_o/L_o > (1/7): # Equation II-1-71 from the Coastal Engineering Manual
        H_o = H_o*(95/700)
    
    # Convert theta_o from radians to degrees
    theta_o = theta_o * 360/(2*pi)
    
    # If the wave angle is impossible, consider shore sheltered from this storm
    if isnan(theta_o):
        theta_o = 0
        H_o = 0
        T_o = 10
    
    return (H_o, T_o, theta_o, L_o, C_o)

def deep2shallow(H_o, T_o, theta_o, L_o, C_o, h_s=20):
    '''Calculates wave shoaling and refraction experienced by offshore, deep-
    water, waves migrating into the nearshore region and returns the resulting
    nearshore wave parameters at a given depth.
    
    Inputs:
        REQUIRED \n
        H_o = Initial wave height (m) \n
        T_o = Initial wave period (sec) \n
        theta_o = Initial wave angle (deg) \n
        L_o = Offshore wave length (m) \n
        C_o = Offshore wave celerity (m/s) \n
        h_s = Nearshore water depth (m) where new waves are calculated \n
    
    Outputs:
        H_s = Nearshore wave height (m) \n
        T_s = Nearshore wave period (sec) \n
        theta_s = Neashore wave angle (deg) \n
        L_s = Nearshore wave length (m) \n
        C_s = Nearshore wave celerity (m/s) \n
    '''

    # Required Constants
    from numpy import pi, tanh, arcsin, sin, cos
    g = 9.8066 #(m/s/s) - Gravitational acceleration
    
    # Convert input angle to radians
    theta_o = theta_o*(2*pi)/360
    
    # Calculate initial wave length, assuming deep water conditions 
    # (ie d/L>0.5)   
    L_s_deep = g*T_o**2/(2*pi)       # Wave length (m)
    
    # Assuming shallow water conditions
    L_s_shallow = (T_o*(g*h_s)**0.5) # Wave length (m)
    
    if h_s/L_s_shallow < 0.05:
        # Shallow water conditions
        L_s = L_s_shallow
        C_s = (g*h_s)**0.5
        
    elif h_s/L_s_deep > 0.5:
        # Deep water conditions
        L_s = L_s_deep
        C_s = ((g*L_o)/(2*pi))**0.5

    else:
        # Intermediate depth wave conditions...this is more complicated
        from scipy.optimize import brentq
        
        # Celerity and wave length for intermediate wave conditions
        # Source: Robert M. Sorensen, "Basic Wave Mechanics for Coastal and 
        # Ocean Engineers" (John Wiley & Sons, 1993), Chapter 2.
        def find_L_i(L, T, h):
            # Equation for the length of an intermediate depth wave needs to be
            # solved iteratively
            zero = L - (g*T**2)/(2*pi) * tanh(2*pi*h/L)
            return zero
        L_s = brentq(find_L_i, 10, max(L_s_deep, L_s_shallow), args=(T_o, h_s))
        C_s = ((g*L_s)/(2*pi) * tanh(2*pi*h_s/L_s))**0.5

    # Converting from deepwater equivalent wave conditions to nearshore wave
    # conditions.
    # Calculate shoaling and refraction coefficients following the methodology
    # presented in 'Water Wave Mechanics for Engineers and Scientists' by
    # Dean & Dalrymple (c) 2000
    K_s = (C_o/C_s)**0.5 # Shoaling coefficient. Eq. 4.116
    
    # Find the offshore equivalent wave angle
    theta_s = arcsin(sin(theta_o)*(C_s/C_o)) # Snell's law, Eq. 4.109
    K_r = (cos(theta_o)/cos(theta_s))**0.5 # Refraction coefficient. Eq. 4.118
    
    # Apply shoaling and refraction coefficients to find new wave height
    H_s = H_o*K_s*K_r # Eq. 4.117
    
    # Wave period is assumed constant
    T_s = T_o
    
    # Convert theta_o from radians to degrees
    theta_s = theta_s * 360/(2*pi) 
    
    return (H_s, T_s, theta_s, L_s, C_s)

def shallow2breaking(H_s, T_s, theta_s, L_s, C_s, h_s, m, A,
                     gamma=0.78, tolerance=1, limit=100):
    '''Returns breaking wave height and depth based on an initial nearshore
    wave condition defined by a wave height, period, and angle originating at a
    given distance from the shoreline. Based on linear wave theory and the 
    assumption of shore-parallel depth contours.
    
    Inputs:
        H_s = Nearshore wave height (m) \n
        T_s = Nearshore wave period (T) \n
        theta_s = Nearshore wave angle from shore perpendicular (deg) \n
        L_s = Nearshore wave length (m) \n
        C_s = Nearshore wave celerity (m/s) \n          
        h_s = Nearshore depth where waves originate (m) \n
        m   = Linear beach-face slope \n
        A   = 2.25*((w**2)/g)**(1/3), from Eq. 15 of Kriebel & Dean (1993). \n
              A parameter governing profile steepness, valid for sand where \n
              0.1mm < d_50 < 0.4mm \n
        gamma = Breaker index, usually taken to be 0.78-1.0. \n
        tolerance = resolution of distance x from shore required for breaking \n
        limit = maximum iterations allowed \n
             
    Outputs:
        H_b = Breaking wave height (m) \n
        h_b = Breaking wave water depth (m) \n
        x_b = Breaker distance offshore (m)
    '''
    x_s = cross_shore_distance(h_s, m, A)
    h_1 = 0.75*h_s
    x_1 = cross_shore_distance(h_1, m, A)
    
    (H_1, T_1, theta_1, L_1, C_1) = deep2shallow(H_s, T_s, theta_s,
                                                 L_s, C_s, h_s=h_1)
    
    def check_breaking(H_1, h_1, x_1, h_s, x_s,
                       gamma=0.78, tolerance=tolerance):  
        if H_1/h_1 > gamma:
            if abs(x_s-x_1) < tolerance:
                h_1 = (h_s+h_1)/2
            else:
                h_1 = 1.175*h_1
        else:
            h_s = h_1
            h_1 = 0.75*h_s
        return (h_s, h_1)
    
    iterations = 0
    while iterations<limit:
        (h_s, h_1) = check_breaking(H_1, h_1, x_1, h_s, x_s)
        x_s = cross_shore_distance(h_s, m, A)
        x_1 = cross_shore_distance(h_1, m, A)
        if (x_s-x_1)<tolerance:
            (H_b, T_b, theta_b, L_b, C_b) = deep2shallow(H_s, T_s, theta_s,
                                                         L_s, C_s, h_s=h_1)
            h_b = h_1
            x_b = (x_s+x_1)/2
            iterations += limit
        (H_s, T_s, theta_s, L_s, C_s) = deep2shallow(H_s, T_s, theta_s,
                                                     L_s, C_s, h_s=h_s)
        (H_1, T_1, theta_1, L_1, C_1) = deep2shallow(H_s, T_s, theta_s,
                                                     L_s, C_s, h_s=h_1)
        iterations += 1
    return (H_b, h_b, x_b) #, T_b, theta_b, L_b, C_b)
    
'''
    h_T = (4/9)*(A**3/m**2) # Depth at which the linear slope is tangent to the
                              concave profile
    x_0 = h_T/(3*m)
'''
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


    