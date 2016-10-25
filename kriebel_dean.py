# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 15:53:51 2016

@author: jdorvinen
"""
# <codecell>
import numpy as np

# <codecell>
# Constants
g = 9.8066 # gravitational acceleration (m/s/s)

# Sediment data
#d_50 = 0.3 # mass-mean sediment grain-size diameter (mm)
w_cm = 2.23 # sediment fall velocity (cm/sec)
w = w_cm/100 # m/sec

# Profile data, based on equilibrium profile of the form 'x=(h/A)**(3/2)'
# where h = the water depth at a distance x offshore from the still-water level
A = 2.25*((w**2)/g)**(1/3) # Eq. 15 'parameter governs profile steepness'
                           # valid for sand where 0.1mm < d_50 < 0.4mm
B = 3.4 #5.5  # Berm height above mean sea-level (meters)
D = 0    # Dune height
#W =  # Width of the back-shore
m = 0.12 #0.10 # Linear beach face slope (m/m)

# Storm data
S = 1.83 #2.74        # given water-level rise ('storm-surge') (meters)
T_d = 66 #8         # Storm duration (hours)
H_b = 6.1*0.78 #4.9*0.78  # Breaking wave height (meters)
gamma = 0.78    # Breaker index, usually taken to be 0.78-1.0.
h_b = H_b/gamma # Breaking depth, assumed to remain constant (meters)

# Active profile width 'x_b', x_0 = the distance from the still-water shoreline
# to the virtual origin of the concave equilibrium profile form, given by x_0 =
# h_T/3m, where h_T is teh depth at which the linear slope is tangent to the 
# concave profile, which may be shown to equal 4A**3/9m**2. Here it is assumed
# that x_0 ~= 0 for simplicity.
x_0 = 0
x_b = x_0+(h_b/A)**(3/2) # Eq. 23

# <codecell>
# Calculate maximum erosion potential 'R_inf' and maximum potential volume
# eroded 'V_inf' based on an equilibrium profile with a linear beach slope.
R_inf = S*(x_b-(h_b/m)) / (B+h_b-(S/2))                        # Eq. 22
V_inf = R_inf*B + (S**2)/(2*m) - (2/5)*(S**(5/2))/(A**(3/2)) # Eq. 24

# Calculate maximum erosion potential 'R_inf' and maximum potential volume
# eroded 'V_inf' based on an equilibrium profile with a dune.

# Dune with no back-shore.
#R_inf = S*(x_b-(h_b/m))/(B+D+h_b-(S/2)) # Eq. 25

# Dune with a wide back-shore.
# R_inf = (S*(x_b-(h_b/m)) - (W*(B+h_b-(S/2)))) / (B+D+h_b-(S/2)) # Eq. 26

# Volume eroded above original sea level in both cases 'V_inf'
#V_inf = R_inf*D + (R_inf+W)*B + (S**2)/(2*m) - (2/5)*(S**(5/2))/(A**(3/2))

# <codecell>
# Time scale of profile response
C_1 = 320 # Empirical coefficient from Kriebel and Dean 1993
# Time scale parameter
T_s = C_1 * \
      ((H_b**(3/2))/(g**(1/2) * A**3)) * \
      (1+(h_b/B)+(m*x_b)/h_b)**(-1) # Eq. 31

# <codecell>
#Beach response to idealized storm surge
alpha = 1/T_s
sigma = np.pi/T_d
beta = 2*sigma/alpha # 2*np.pi*(T_s/T_d)

# Eq. 10
# R_t/R_inf = 0.5*(1 - \
#                  (beta**2/(1+beta**2))*np.exp(-(2*sigma*t)/beta) - \
#                  (1/(1+beta**2))*(np.cos(2*sigma*t)+beta*np.sin(2*sigma*t)) \
#                 )

# Setting time derivative of Eq. 10 to zero leads to Eq. 12, where t_max is
# the time at which maximum erosion will take place.
import scipy.optimize as opt

def func(t_max):
    zero = np.cos(2*sigma*t_max) - \
           (1/beta)*np.sin(2*sigma*t_max) - \
           np.exp(-(2*sigma*t_max)/beta)
    return zero

# <codecell>
t_max = opt.brentq(func,
                   a=T_d/2,
                   b=T_d)

R_max = R_inf*0.5*(1-np.cos(2*t_max))

#print(R)