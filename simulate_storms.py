# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# In[]
import os
import numpy as np
import scipy.stats as st
import pandas as pd
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
    '''Maximum Hsig is modeled by a generalized pareto distribution'''
    c = -0.17228
    loc = 3.1125
    scale = 1.2256
    #rnv_copulaized = hsig_storm_len_copula(rnv, storm_len) 
    hsig = int(round(loc + ((1-rnv)**-c - 1) * scale/c))   
    #hsig = st.genpareto.rvs(c=c, loc=loc, scale=scale, discrete=True)
    return hsig

def get_a_hsig(rnv, hsig):
    '''Average Hsig is modeled by a generalized pareto distribution'''
    c = -0.21576
    loc = 3.0766
    scale = 0.59362
    rnv_copulaized = a_hsig_hsig_copula(rnv, hsig) 
    a_hsig = int(round(loc + ((1-rnv_copulaized)**-c - 1) * scale/c)) 
    #hsig = st.genpareto.rvs(c=c, loc=loc, scale=scale, discrete=True)
    return a_hsig

def get_tps(rnv, hsig):
    '''Modeled as a Frechet distribution'''
    c = 10.503
    loc = -4.9823
    scale = 13.506
    #rnv_copulaized = tps__hsig_copula(rnv, hsig)
    tps = loc - scale*-(abs(np.log(rnv))**(-1/c))
    #tps2 = st.genextreme.rvs(c=-0.08, loc=8.58, scale=1.31, discrete=False)
    #tps = [tps1, tps2]
    return tps

def get_a_tps(rnv, tps):
    '''Modeled as a Frechet distribution'''
    c = 12.409
    loc = -5.2135
    scale = 13.723
    #rnv_copulaized = tps__hsig_copula(rnv, hsig)
    a_tps = loc - scale*-(abs(np.log(rnv))**(-1/c))    
    #a_tps2 = st.genextreme.rvs(c=-0.08, loc=8.55, scale=1.12, discrete=False)
    #a_tps = [a_tps1, a_tps2]
    return a_tps
    
def get_tide(rnv, storm_len):
    '''Maximum Hsig is modeled by a generalized pareto distribution'''
    c = -0.17228
    loc = 3.1125
    scale = 1.2256
    #rnv_copulaized = hsig_storm_len_copula(rnv, storm_len) 
    tide = 1.0 / (1 + exp(-1.0 / b * (np.random.normal() - a)))
    #hsig = st.genpareto.rvs(c=c, loc=loc, scale=scale, discrete=True)
    return hsig
    
SWEL
johnsonsb(a=-0.95, b=1.41, loc=-1.12, scale=1.90)<
KstestResult(statistic=0.029574404259407161, pvalue=0.9901000531575973)<br>

\gamma = -1.1319
\rho = 1.5523
\lambda = 2.098
\epsilon = -1.2726

# In[]
i = 0
storm = []
while i < 100000:
    strm = get_a_tps(np.random.random(),5)
    storm.append(strm)
    i += 1

tps1 = [storm[i][0] for i in range(100000)]
tps2 = [storm[i][1] for i in range(100000)]

tp1 = pd.Series(data=tps1)
tp2 = pd.Series(data=tps2)
    
# In[]
i = 0
Time = [0]
max_time = 100*8766

 
johnsonsb = johnsonsb_gen(a=0.0, b=1.0, name='johnsonsb')

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
