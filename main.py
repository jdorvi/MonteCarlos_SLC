# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 14:31:54 2016

@author: jdorvinen
"""
from wave_breaking import input2deep, deep2shallow, shallow2breaking, wave_setup
from kriebel_dean import kriebel_dean, recovery

# Define input parameters
w_cm = 10 # Sediment fall velocity (cm/s)
B = 2 #berm height (meters)
D = 1 #Dune height (meters)
W = 10 #width of the back-shore (meters)
m = 0.05 #Linear beach face slope (m/m)
A = 0.175 # 2.25*((w**2)/g)**(1/3), from Eq. 15 of Kriebel & Dean (1993).
        # A parameter governing profile steepness, valid for sand where
        # 0.1mm < d_50 < 0.4mm
SLR_rate = 0.0035 # sea level rise rate (m/year)
def main():
    ''' Bringing it all together'''
    import pandas as pd
    
    V_eroded = 0 # Initially eroded volume
    i = 1
    storms = pd.read_pickle('temp_storms.npy')
    
    volume = [V_eroded]
    time = [0]
    for storm in storms.T:
        V_max = V_eroded
        # Storm parameters
        H_i = storms['hsig'][storm]
        T_i = storms['hsig'][storm]
        theta_i = 0 #storms['hsig'][storm]
        # Tidal anomoly
        tide = storms['tide'][storm]
        # Storm length and interim
        interim = storms['interim'][storm]
        length = storms['length'][storm]

        # Input storm conditions to offshore storm conditions
        (H_o, T_o, theta_o, L_o, C_o)=input2deep(H_i, T_i, theta_i, 52.4)
        
        # Bring offshore conditions to the nearshore
        (H_s, T_s, theta_s, L_s, C_s)=deep2shallow(H_o, T_o, theta_o, L_o, C_o)

        # Find breaking wave characteristics
        (H_b, h_b, x_b)=shallow2breaking(H_s, T_s, theta_s, L_s, C_s, h_s=20, m=m, A=A)
        
        # Calculate wave setup
        setup = wave_setup(H_b)
        
        # Calculate surge
        surge = setup+tide
        
        # Add SLR
        if i > 2:
            surge += SLR_rate*time[i-2]/(24*365.25)
        
        # Calculate recovery
        V_recovered = recovery(V_max, interim)
        
        # Calculate erosion
        (V_eroded, R_eroded)=kriebel_dean(w_cm, B, D, W, m, surge, length, H_b)
        if V_eroded < V_recovered:
            V_eroded = V_recovered
            
        # add to lists
        volume.append(V_recovered)
        volume.append(V_eroded)
        time.append(time[i-1]+interim)
        time.append(time[i]+length)
        
        # Increment time index
        i += 2
        
    return (time, volume)
    
if __name__ == '__main__':
    (time, volume) = main()
