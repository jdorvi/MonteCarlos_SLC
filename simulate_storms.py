# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# In[]
import os
import numpy as np
import scipy.stats as st
# In[]
outpath = "C:/Users/jdorvinen/Documents/Jared/Projects/East Hampton/"
outfile = "test1.csv"
outputfile = os.path.join(outpath, outfile)
# In[]
def get_interim(rnv):
    '''Estimate time between storms based on a folded Cauchy distribution'''
    c = 0.63
    loc = 25.94
    scale = 137.16
    interim = int(round(abs(scale * np.tan(np.pi * (rnv-1/2)) + loc)))
    #interim2 = st.foldcauchy.rvs(c=c, loc=loc, scale=scale, discrete=True)
    #interim = [interim1, interim2]
    return interim

def get_storm_len(rnv):
    '''Storm length is modeled by a generalized pareto distribution'''
    c = -0.054     # -0.02
    loc = 1.9853   #  3.0
    scale = 15.852 # 14.16 
    storm_len = int(round(loc + (((1-rnv)**-c)-1)*scale/c)) 
    #storm_len2 = st.genpareto.rvs(c=c, loc=loc, scale=scale, discrete=True)
    #storm_len = [storm_len1, storm_len2]
    return storm_len

def get_hsig(rnv, storm_len):
    c = -0.17228
    loc = 3.1125
    scale = 1.2256
    hsig = int(round(loc + ((1-rnv)**-c - 1) * scale/c))   
    #hsig = st.genpareto.rvs(c=c, loc=loc, scale=scale, discrete=True)
    return hsig

def get_a_hsig(rnv, hsig):
    a = 0.22
    b = 2.34
    c = 7.13
    loc = 3.03
    scale = 1.09
    rnv_copulaized = a_hsig__hsig_copula(rnv, hsig) 
    a_hsig = int(round(loc + ((1-rnv_copulaized)**-c - 1) * scale/c)) 
    #a_hsig = st.genexpon.rvs(a=a, b=b, c=c, loc=loc, scale=scale, discrete=False)
    return a_hsig

def get_tps(rnv, hsig):
    '''Modeled as a Frechet right distribution'''
    c = 1.5975
    loc = 6.5955
    scale = 3.1261
    rnv_copulaized = tps__hsig_copula(rnv, hsig)
    tps = loc + scale*np.log(rnv)**(-1/c)
    #tps = st.weibull_min.rvs(mu=mu, loc=loc, scale=scale, discrete=False)
    return tps

def get_a_tps(rnv, tps):
    mu = 0.22
    loc = 5.90
    scale = 0.61
    
    
    #a_tps = st.recipinvgauss.rvs(mu=0.22, loc=5.90, scale=0.61, discrete=False)
    #a_tps = tps-rnv*3
    return a_tps
    
# In[]
i = 0
Time = [0]
max_time = 100*8766
# In[]
%%timeit
i = 0
Time = [0]
with open(outputfile, 'w') as outfile:
    while Time[i] < max_time:
        #rnv = [np.random.random() for variable in range(6)]
        interim = get_interim(rnv[0], Time[i])
        storm_len = get_storm_len(rnv[1])
        hsig = get_hsig(rnv[2], storm_len)
        a_hsig = get_a_hsig(rnv[3], hsig)
        tps = get_tps(rnv[4], hsig)
        a_tps = get_a_tps(rnv[5], tps)
        
        raw_line = "{:>8},{:>8},{:>8.2f},{:>8.2f},{:>8.2f},{:>8.2f}\n"
        line = raw_line.format(interim, storm_len, hsig, a_hsig, tps, a_tps)
        outfile.write(line)
        
        Time.append(Time[i]+interim+storm_len)
        i += 1


# In[]
import scipy as sp
import scipy.stats as st
st.foldcauchy.rvs()
