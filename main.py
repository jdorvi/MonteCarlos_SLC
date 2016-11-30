# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 14:31:54 2016

@author: jdorvinen
"""
from wave_breaking import input2deep, deep2shallow, shallow2breaking, wave_setup
from kriebel_dean import kriebel_dean, recovery

# Define input parameters
w_cm = 17.42 # Sediment fall velocity (cm/s)
B = 2 # berm height (meters)
D = 0 # Dune height (meters)
W = 0 # width of the back-shore (meters)
m = 0.1 # Linear beach face slope (m/m)
A = 0.328 # 2.25*((w**2)/g)**(1/3), from Eq. 15 of Kriebel & Dean (1993).
        # A parameter governing profile steepness, valid for sand where
        # 0.1mm < d_50 < 0.4mm
SLR_rate = 0 # sea level rise rate (m/year) (this really doesn't work well)
MAX_INTER = 24 # one half the maximum time between points during recovery (hrs)
T_a = 400 # Recovery characteristic time (hrs)
STORMS = 'temp_storms3.npy' # Simulated storms
SHOREANGLE = -155 # Shoreline angle relative to a north facing shoreline (degs)
BUOY_DEPTH = 52.4 # Water depth at offshore wave buoy (m)
NEARSHORE_DEPTH = 20 # Water depth nearshore where wave breaking calcs begin(m)

def main():
    ''' Bringing it all together'''
    import pandas as pd
    
    V_eroded = 0 # Initially eroded volume
    R_eroded = 0 # Initially eroded distance
    i = 1 # Start timestep index at one
    storms = pd.read_pickle(STORMS) # Read in saved storms
    
    # Initialize volume, distance, and time lists which will contain the 
    # outputs of this program.
    volume = [V_eroded]
    distance = [R_eroded]
    time = [0]

    # Define function to break up interval between storm events and create
    # smoother plots
    def inter(num):
        """ Calculates number of intervals of length MAX_INTER and remainder
        for a storm interval of a given length, num """
        if num/(2*MAX_INTER) > 1:
            intervals = int(((num/MAX_INTER)//2)*2)
            remainder = num%(intervals*MAX_INTER)
        else:
            intervals = 0
            remainder = num
        return (intervals, remainder)

    # Walk through all saved storms and simulate their corresponding erosion
    for storm in storms.T:
        V_max = V_eroded
        R_max = R_eroded
        
        # Storm parameters
        H_i = storms['hsig'][storm] # Significant wave height
        T_i = storms['tps'][storm] # Dominant wave period
        
        # Wave direction relative to shore-normal
        theta_i = storms['direction'][storm] - SHOREANGLE

        # Tidal anomoly
        tide = storms['tide'][storm]

        # Storm length and interim
        interim = storms['interim'][storm]
        length = storms['length'][storm]

        # Calculate subintervals for recovery
        (intervals, remainder) = inter(interim)

        # Input storm conditions to offshore storm conditions
        (H_o, T_o, theta_o, L_o, C_o)=input2deep(H_i, T_i, theta_i, BUOY_DEPTH)
        
        # Bring offshore conditions to the nearshore
        (H_s, T_s, theta_s, L_s, C_s)=deep2shallow(H_o, T_o, theta_i, L_o, C_o)

        # Find breaking wave characteristics
        (H_b, h_b, x_b)=shallow2breaking(H_s,
                                         T_s,
                                         theta_s,
                                         L_s,
                                         C_s,
                                         h_s=NEARSHORE_DEPTH,
                                         m=m,
                                         A=A)
        
        # Calculate wave setup
        setup = wave_setup(H_b)
        
        # Calculate surge
        surge = setup+tide
        
        # Calculate recovery before storm for each interval
        for j in range(intervals):
            V_recovered = recovery(V_max, (1+j)*MAX_INTER, T_a=T_a)
            R_recovered = recovery(R_max, (1+j)*MAX_INTER, T_a=T_a)
            volume.append(V_recovered)
            distance.append(R_recovered)        
            time.append(time[i-1]+MAX_INTER)
            i += 1
        
        # Calculate final recovery before storm
        V_recovered = recovery(V_max, interim, T_a=T_a)
        R_recovered = recovery(R_max, interim, T_a=T_a)
        volume.append(V_recovered)
        distance.append(R_recovered)        
        time.append(time[i-1]+remainder)
        i += 1
        
        # Calculate hypothetical recovery during storm
        V_recovered = recovery(V_recovered, length, T_a=T_a)
        R_recovered = recovery(R_recovered, length, T_a=T_a)
        
        # Add SLR
        if i > 2:
            surge += (SLR_rate*time[i-1]/(24*365.25))
            
        # Calculate hypothetical erosion during storm
        (V_eroded, R_eroded, V_inf, R_inf)=kriebel_dean(w_cm,
                                                        B,
                                                        D,
                                                        W,
                                                        m,
                                                        surge,
                                                        length,
                                                        H_b)
        
        # Check if storm erosion is greater than current erosion, choose the 
        # larger of the two
        '''
        if V_eroded < V_recovered:
            V_eroded = V_recovered
        if R_eroded < R_recovered:
            R_eroded = R_recovered
        '''
        # Sum current and storm erosion 
        if V_eroded > 0:
            V_eroded += V_recovered
            R_eroded += R_recovered
        else:
            V_eroded = V_recovered
            R_eroded = R_recovered
            
        # Add eroded volume and distance to respective lists
        volume.append(V_eroded)
        distance.append(R_eroded)
        time.append(time[i-1]+length)
        
        # Increment time index
        i += 1
        
    return (time, volume, distance)
    
if __name__ == '__main__':
    (time, volume, distance) = main()
