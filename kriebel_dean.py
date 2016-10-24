# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 15:53:51 2016

@author: jdorvinen
"""

# Constants
g = 9.81 # gravitational acceleration (m/s/s)

# Sediment data
d_50 = # mass-mean sediment grain-size diameter
w = # sediment fall velocity

# Profile data, based on equilibrium profile of the form 'x=(h/A)**(3/2)'
# where h = the water depth at a distance x offshore from the still-water level
A = 2.25*((w**2)/g)**(1/3) # Eq. 15 'parameter governs profile steepness'
                           # valid for sand where 0.1mm < d_50 < 0.4mm
B = # Berm height above mean sea-level
m = # Linear beach face slope 

# Storm data
S = # given water-level rise ('storm-surge')
T_d = # 
H_b = # Breaking wave height
gamma = 0.78 # Breaker index, usually taken to be 0.78-1.0.
h_b = H_b/gamma # Breaking depth, assumed to remain constant. 

# Active profile width 'x_b', x_0 = the distance from the still-water shoreline
# to the virtual origin of the concave equilibrium profile form, given by x_0 =
# h_T/3m, where h_T is teh depth at which the linear slope is tangent to the 
# concave profile, which may be shown to equal 4A**3/9m**2. Here it is assumed
# that x_0 ~= 0 for simplicity.
x_0 = 0
x_b = x_0+(h_b/A)**(3/2) # Eq. 23

# Calculate maximum erosion potential 'R_inf' and maximum potential volume
# eroded 'V_inf' based on an equilibrium profile with a linear beach slope.
R_inf = S*(x_b-(h_b/m))/(B+h_b+(S/2)) # Eq. 22 
V_inf = R_inf*B + (S**2)/(2*m) - (2/5)*(S**(5/2))/(A**(3/2)) # Eq. 24


# Calculate maximum erosion potential 'R_inf' and maximum potential volume
# eroded 'V_inf' based on an equilibrium profile with a dune.
D =  # Dune height
W =  # Width of the back-shore

# Dune with no back-shore.
R_inf = S*(x_b-(h_b/m))/(B+D+h_b-(S/2)) # Eq. 25

# Dune with a wide back-shore.
R_inf = (S*(x_b-(h_b/m)) - (W*(B+h_b-(S/2)))) / (B+D+h_b-(S/2)) # Eq. 26

# Volume eroded above original sea level in both cases 'V_inf'
V_inf = R_inf*D + (R_inf+W)*B + (S**2)/(2*m) - (2/5)*(S**(5/2))/(A**(3/2))