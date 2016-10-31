# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 14:31:54 2016

@author: jdorvinen
"""
from wave_breaking import input2deep, deep2shallow, shallow2breaking, wave_setup
from kriebel_dean import kriebel_dean, recovery

# Define input parameters
w_cm = 4.86 # Sediment fall velocity (cm/s)
B = 10 # berm height (meters)
D = 0 # Dune height (meters)
W = 0 # width of the back-shore (meters)
m = 0.07 # Linear beach face slope (m/m)
A = 0.14 # 2.25*((w**2)/g)**(1/3), from Eq. 15 of Kriebel & Dean (1993).
        # A parameter governing profile steepness, valid for sand where
        # 0.1mm < d_50 < 0.4mm
SLR_rate = 0 # sea level rise rate (m/year)
def main():
    ''' Bringing it all together'''
    import pandas as pd
    
    V_eroded = 0 # Initially eroded volume
    R_eroded = 0
    i = 1
    storms = pd.read_pickle('temp_storms2.npy')
    
    volume = [V_eroded]
    distance = [R_eroded]
    time = [0]

    for storm in storms.T:
        V_max = V_eroded
        R_max = R_eroded
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
        V_recovered = recovery(V_max, interim, T_a=400)
        R_recovered = recovery(R_max, interim, T_a=400)
        # Calculate erosion
        (V_eroded, R_eroded, V_inf, R_inf)=kriebel_dean(w_cm, B, D, W, m, surge, length, H_b)
        #print(V_eroded)
        if V_eroded < V_recovered:
            V_eroded = V_recovered
        if R_eroded < R_recovered:
            R_eroded = R_recovered
            
        # add to lists
        volume.append(V_recovered)
        volume.append(V_eroded)
        distance.append(R_recovered)
        distance.append(R_eroded)
        time.append(time[i-1]+interim)
        time.append(time[i]+length)
        
        # Increment time index
        i += 2
        
    return (time, volume, distance)
    
if __name__ == '__main__':
    (time, volume, distance) = main()
