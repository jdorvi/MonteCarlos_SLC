# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 11:54:17 2016

@author: jdorvinen
"""
def fall_velocity(d_50, rho_s=2650, rho_w=1028, a=1, b=1, c=0.9):
    ''' Calculate sediment fall velocity based the method presented in, 
    
    Wu, W., and Wang, S.S.Y., 'Formulas for sediment porosity and settling
    velocity' J. Hydraul. Eng., 2006, 132(8): 858-862
     
    Inputs:
        REQUIRED
        d_50 = Sediment mass mean diameter grainsize (mm)
        OPTIONAL
        rho_s = 2650 (kg/m**3), sediment density (quartz value as default)
        rho_w = 1028 (kg/m**3), sea water density (Gill's (1982) reference as
                default. Where salinity = 35 PPT, temp = 5 degC, and pressure 
                = 0 dbar)
        a,b,c = Lengths of the longest, intermediate, and shortest axes of the 
                sediment particle. (dimensions not important)
    
    Returns:
        w = sediment fall velocity (cm/s) '''

    from numpy import exp
    
    d = d_50 / 1000 # sediment mass mean diameter grainsize (m)
    g = 9.8066      # gravitational acceleration (m/s/s)
    s = rho_s/rho_w # relative density
    dv = 0.0016193  # dynamic viscosity (kg/(m*s)), S=35PPT, temp=5degC, P=1atm
    kv = dv/rho_w   # m**2/s
    
    S_f = c/((a*b)**0.5) # Corey shape factor 
    M = 53.5*exp(-0.65*S_f)
    N = 5.65*exp(-2.5*S_f)
    n = 0.7+0.9*S_f
    D_s = d*((rho_s/(rho_w-1))*(g/kv**2))**(1/3)
    
    # Calculate fall-velocity
    # Cut-off of 0.2557 mm was chosen based on approximate equivalence of 
    # methods at this grain-sized diameter
    if d_50 > 0.2557:
        '''Wu et al. 2006'''
        w = (M*kv)/(N*d) * \
            (((0.25+(4*N*D_s**3)/(3*M**2))**(1/n))**0.5 - 0.5)**n
    elif d_50 <= 0.2557:
        '''Stokes equation'''
        w = 1/18 * (s-1)*(g*d**2)/kv

    w *= 100 # convert m/s to cm/s
    
    return w
    