# -*- coding: utf-8 -*-
"""
This is a script for simulating a series of random storm events from 
distributions fitted to historical observations. 

author: j. dorvinen
email: jdorvi_at_gmail_dot_com
date: 10/13/2016
"""

# <codecell>
import os
import numpy as np
import scipy.stats as st
import pandas as pd
# <codecell>
outpath = "/home/rannikko/git/"
outfile = "test1.csv"
outputfile = os.path.join(outpath, outfile)
# <codecell>
def get_interim(rnv):
    '''Estimate time between storms based on a folded Cauchy distribution. The
    given configuration gives a period of hours as an integer.'''
    c = 0.63
    loc = 25.94
    scale = 137.16
    interim = int(round(abs(scale * np.tan(np.pi * (rnv-1/2)) + loc)))
    #interim2 = st.foldcauchy.rvs(c=c, loc=loc, scale=scale, discrete=True)
    #interim = [interim1, interim2]
    return interim

def get_storm_len(rnv):
    '''Storm length is modeled by a generalized Pareto distribution. The given
    configuration produces a storm of some hours length as an integer.'''
    # distribution fitting paramters
    c = -0.054     # -0.02
    loc = 1.9853   #  3.0
    scale = 15.852 # 14.16
    # storm_len is sampled via the fitted quantile function.
    storm_len = int(round(loc + (((1-rnv)**-c)-1)*scale/c)) 
    #storm_len2 = st.genpareto.rvs(c=c, loc=loc, scale=scale, discrete=True)
    #storm_len = [storm_len1, storm_len2]
    return storm_len

def get_hsig(rnv):
    '''Maximum Hsig is modeled by a generalized Pareto distribution. In the 
    current configuration Hsig is measured in meters.'''
    c = -0.17228
    loc = 3.1125
    scale = 1.2256
    #rnv_copulaized = hsig_storm_len_copula(rnv, storm_len) 
    hsig = loc + ((1-rnv)**-c - 1) * scale/c  
    #hsig = st.genpareto.rvs(c=c, loc=loc, scale=scale, discrete=True)
    return hsig

def get_a_hsig(rnv):
    '''Average Hsig is modeled by a generalized Pareto distribution'''
    c = -0.21576
    loc = 3.0766
    scale = 0.59362
    #rnv_copulaized = a_hsig_hsig_copula(rnv, hsig) 
    a_hsig = loc + ((1-rnv)**-c - 1) * scale/c 
    #hsig = st.genpareto.rvs(c=c, loc=loc, scale=scale, discrete=True)
    return a_hsig

def get_tps(rnv):
    '''Modeled as a Frechet distribution'''
    c = 10.503
    loc = -4.9823
    scale = 13.506
    #rnv_copulaized = tps__hsig_copula(rnv, hsig)
    tps = loc - scale*-(abs(np.log(rnv))**(-1/c))
    #tps2 = st.genextreme.rvs(c=-0.08, loc=8.58, scale=1.31, discrete=False)
    #tps = [tps1, tps2]
    return tps

def get_a_tps(rnv):
    '''Modeled as a Frechet distribution'''
    c = 12.409
    loc = -5.2135
    scale = 13.723
    #rnv_copulaized = tps__hsig_copula(rnv, hsig)
    a_tps = loc - scale*-(abs(np.log(rnv))**(-1/c))    
    #a_tps2 = st.genextreme.rvs(c=-0.08, loc=8.55, scale=1.12, discrete=False)
    #a_tps = [a_tps1, a_tps2]
    return a_tps
    
def get_tide(rnv):
    '''Maximum tide is modeled by as an independent random variable following
    a JohnsonSB distribution'''
    a = -1.1319   #-0.95
    b = 1.5523    # 1.41
    loc = -1.2726 #-1.12
    scale = 2.098 # 1.90
    tide = scale*1.0 / (1 + np.exp(-1.0 / b * (np.random.normal() - a)))+loc
    #tide2 = st.johnsonsb.rvs(a=a, b=b, loc=loc, scale=scale, discrete=False)
    #tide = [tide1, tide2]
    return tide

# <codecell>
%%time
import chaospy as cp

n_samples = 1000000

joint = cp.J(cp.Uniform(lo=0,up=1),
             cp.Uniform(lo=0,up=1),
             cp.Uniform(lo=0,up=1),
             cp.Uniform(lo=0,up=1),
             cp.Uniform(lo=0,up=1)
            )
clayton = cp.Clayton(joint, theta=2) #why 2?
samples = clayton.sample(size=n_samples)
independent = [cp.Uniform(lo=0,up=1).sample(size=n_samples),
               cp.Uniform(lo=0,up=1).sample(size=n_samples)
              ]
rnv = {'length':samples[0],
       'hsig':samples[1],
       'a_hsig':samples[2],
       'tps':samples[3],
       'a_tps':samples[4],
       'interim':independent[0], #needs to be reformulated
       'tide':independent[1]
      }

storms = {'length':[get_storm_len(rnv['length'][i]) for i in range(n_samples)],
          'hsig':[get_hsig(rnv['hsig'][i]) for i in range(n_samples)],
          'a_hsig':[get_a_hsig(rnv['a_hsig'][i]) for i in range(n_samples)],                  
          'tps':[get_tps(rnv['tps'][i]) for i in range(n_samples)],                  
          'a_tps':[get_a_tps(rnv['a_tps'][i]) for i in range(n_samples)],
          'interim':[get_interim(rnv['interim'][i]) for i in range(n_samples)],
          'tide':[get_tide(rnv['tide'][i]) for i in range(n_samples)]
         }
                  
# <codecell>
i = 0
storm = []
while i < 100000:
    strm = get_tide()
    storm.append(strm)
    i += 1

tps1 = [storm[i][0] for i in range(100000)]
tps2 = [storm[i][1] for i in range(100000)]

tp1 = pd.Series(data=tps1)
tp2 = pd.Series(data=tps2)
    
# <codecell>
i = 0
Time = [0]
max_time = 1000*8766

# <codecell>
%%timeit
i = 0
Time = [0]
with open(outputfile, 'w') as outfile:
    header = " interim,  duratn,    tide,    hsig,  a_hsig,     tps,   a_tps\n"
    outfile.write(header)
    while Time[i] < max_time:
        rnv = [np.random.random() for variable in range(6)]
        interim = get_interim(rnv[0])
        storm_len = get_storm_len(rnv[1])
        hsig = get_hsig(rnv[2], storm_len)
        a_hsig = get_a_hsig(rnv[3], hsig)
        tps = get_tps(rnv[4], hsig)
        a_tps = get_a_tps(rnv[5], tps)
        tide = get_tide()
        raw_line = "{:>8},{:>8},{:>8.3f},{:>8.2f},{:>8.2f},{:>8.2f},{:>8.2f}\n"
        line = raw_line.format(interim,
                               storm_len,
                               tide,
                               hsig,
                               a_hsig,
                               tps,
                               a_tps)
        outfile.write(line)
        
        Time.append(Time[i]+interim+storm_len)
        i += 1

# <codecell>
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(1,figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')
cmpair = plt.cm.get_cmap('Paired')
colors = [np.random.random() for i in range(0,n_samples)]
for storm in storms:
    ax.scatter(storms['hsig'][:],
               storms['tps'][:],
                storms['length'][:],
                'o',
                c=colors,
                cmap=cmpair,
                alpha=0.5)
fig.show()