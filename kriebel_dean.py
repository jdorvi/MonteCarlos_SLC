# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 15:53:51 2016

@author: jdorvinen
"""

import numpy as np

def kriebel_dean(w_cm, B, D, W, m, S, T_d, H_b, gamma=0.78):
    '''Calculates storm erosion based on the method presented in,
    
    Kriebel, D.L., and Dean, R.G., 'Convolution method for time-dependent 
    beach-profile response' J. Waterway, Port, Coastal, Ocean Eng., 1993,
    119(2): 204-226

    Inputs:
        REQUIRED
        w_cm = sediment fall velocity (cm/s)
        B    = berm height (meters)
        D    = Dune height (meters)
        W    = width of the back-shore (meters)
        m    = Linear beach face slope (m/m)
        S    = Water-level rise ('storm-surge') (meters)
        T_d  = Storm duration (hours)
        H_b  = Breaking wave height (meters)
        OPTIONAL
        gamma = Breaker index, usually taken to be 0.78-1.0 
    
    Returns:
        R_max = Maximum shoreline erosion (meters)'''
        
    # Constants
    g = 9.8066 # gravitational acceleration (m/s/s)
    
    # Sediment data
    #d_50 = 0.3  # mass-mean sediment grain-size diameter (mm)
    w_cm = w_cm  # sediment fall velocity (cm/s)
    w = w_cm/100 # m/sec
    
    # Profile data
    # Based on equilibrium profile of the form 'x=(h/A)**(3/2)', where h = the
    # water depth at a distance x offshore from the still-water level
    A = 2.25*((w**2)/g)**(1/3) # Eq. 15 'parameter governs profile steepness'
                               # valid for sand where 0.1mm < d_50 < 0.4mm
    B = B   # Berm height above mean sea-level (meters)
    D = D   # Dune height (meters)
    W = W   # Width of the back-shore (meters)
    m = m   # Linear beach face slope (m/m)
    
    # Storm data
    S = S           # given water-level rise ('storm-surge') (meters)
    T_d = T_d       # Storm duration (hours)
    gamma = gamma   # Breaker index, usually taken to be 0.78-1.0.
    H_b = H_b       # Breaking wave height (meters)
    h_b = H_b/gamma # Breaking depth, assumed to remain constant (meters)
    
    # Active profile width 'x_b', x_0 = the distance from the still-water 
    # shoreline to the virtual origin of the concave equilibrium profile form,
    # given by x_0 = h_T/3m, where h_T is teh depth at which the linear slope
    # is tangent to the concave profile, which may be shown to equal 
    # 4A**3/9m**2.
    h_T = (4/9)*(A**3/m**2)  # Eq. 16b_1
    x_0 = h_T/(3*m)          # Eq. 16b_2
    x_b = x_0+(h_b/A)**(3/2) # Eq. 23

    # Calculate erosion potential
    # Maximum erosion potential, 'R_inf', and maximum potential volume eroded,
    # 'V_inf', based on an equilibrium profile with a linear beach slope.
    #R_inf = S*(x_b-(h_b/m)) / (B+h_b-(S/2))                      # Eq. 22
    #V_inf = R_inf*B + (S**2)/(2*m) - (2/5)*(S**(5/2))/(A**(3/2)) # Eq. 24
    
    # Calculate maximum erosion potential 'R_inf' and maximum potential volume
    # eroded 'V_inf' based on an equilibrium profile with a dune.
    
    # Dune with no back-shore.
    # R_inf = S*(x_b-(h_b/m))/(B+D+h_b-(S/2)) # Eq. 25
    
    # Dune with a wide back-shore. 
    R_inf = (S*(x_b-(h_b/m)) - (W*(B+h_b-(S/2)))) / (B+D+h_b-(S/2)) # Eq. 26
    
    # Volume eroded
    # V_inf = R_inf*D + (R_inf+W)*(B-S) # Eq. 27 --> used in K&D examples
    # Volume eroded above original sea level #Eq. 28
    # V_minf = R_inf*D +(R_inf+W)*B+(S**2)/(2*m)-(2/5)*(S**(5/2))/(A**(3/2)) 
    
    # Calculate erosion timescale
    # Time scale of profile response
    C_1 = 320 # Empirical coefficient from Kriebel and Dean 1993
    # Time scale parameter # Eq.31 (sec)
    T_sec = ((H_b**(3/2))/(g**(1/2) * A**3)) / (1+(h_b/B)+(m*x_b)/h_b)
    T_s = C_1*T_sec/3600 # convert seconds to hours
    
    # Combine erosion potential and timescale
    #Beach response to idealized storm surge
    alpha = 1/T_s
    sigma = np.pi/T_d
    beta = 2*sigma/alpha # 2*np.pi*(T_s/T_d)
    
    # Eq. 10
    # R_t/R_inf=0.5*(1 - \
    #               (beta**2/(1+beta**2))*np.exp(-(2*sigma*t)/beta) - \
    #               (1/(1+beta**2))*(np.cos(2*sigma*t)+beta*np.sin(2*sigma*t)))
    
    # Setting time derivative of Eq. 10 to zero leads to Eq. 12, where t_max is
    # the time at which maximum erosion will take place.
    def find_t_max(t_max):
        zero = np.cos(2*sigma*t_max) - \
               (1/beta)*np.sin(2*sigma*t_max) - \
               np.exp(-(2*sigma*t_max)/beta) # Eq. 12
        return zero
    # This can then be solved iteratively to find the time at which maximum
    # erosion occurs, 't_max' (hrs)
    import scipy.optimize as opt
    t_max = opt.brentq(find_t_max,
                       a=T_d/2,
                       b=T_d)
    
    # Finally calculate maximum shoreline recession and volumetric erosion for
    # the given storm parameters.
    R_max = R_inf*0.5*(1-np.cos(2*sigma*t_max)) # Eq. 13
    #V_max = V_inf*(R_max/R_inf)
    
    # Turn this block on if need to debug
    '''
    print("R_max:       {:.1f} (m)".format(R_max))
    print("R_inf:       {:.1f} (m)".format(R_inf))
    print("R_max/R_inf: {:.2f}".format(R_max/R_inf))
    print("V_max:       {:.1f} (m**3/m)".format(V_max))
    print("V_inf:       {:.1f} (m**#/m)".format(V_inf))
    print("T_s:         {:.2f} (h)".format(T_s))
    print("t_max:       {:.1f} (h)".format(t_max))
    print("A:           {:.3f}".format(A))
    print("alpha:       {:.3f} (1/h)".format(alpha))
    print("beta:        {:.3f}".format(beta))
    print("sigma:       {:.3f}".format(sigma))
    '''
    
    return R_max
