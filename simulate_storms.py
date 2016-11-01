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
#import scipy.stats as st
import pandas as pd
import scipy.optimize as opt
# <codecell>
outpath = os.getcwd() #"/home/rannikko/git/"
outfile = "test4.csv"
outputfile = os.path.join(outpath, outfile)
# <codecell>
def interim_func(Gi, te, rnv):
    # Model data fit: alpha=1.498, beta=-0.348, gamma=1.275
    # Callaghan et al. used: alpha=21.46, beta=1.08, gamma=1.07
    a = 17.976 # 12*1.498  # {0}
    b = -4.176 # 12*-0.348 # {1}
    c = 15.3   # 12*1.275  # {2}
    rnv = rnv
    te = te
    w = 2*np.pi
    f = '1 - rnv - \
         np.exp(-({0}*w*Gi + \
                  {1}*(np.cos(w*te) - np.cos(w*(te+Gi))) - \
                  {2}*(np.sin(w*te) - np.sin(w*(te+Gi))))/w)'

    #print(f.format(a,b,c,rnv,te))

    s = eval(f.format(a,b,c))
    return s


def get_interim(rnv,te):
    '''Estimate time between storms based on a non-homogeneous Poisson 
    distribution. The given configuration gives a period of hours as an 
    integer.'''
    try:
        interim_yrs = opt.zeros.brentq(interim_func, 0, 1, 
                                       args=(te,rnv),
                                       xtol=2e-12,
                                       rtol=8.8817841970012523e-16,
                                       maxiter=100,
                                       full_output=False,
                                       disp=True)
        interim = int(round(interim_yrs*365.25*24))
    except ValueError:
        interim = np.nan
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

def get_direction(rnv1, rnv2):
    '''Wave direction is sampled from an empirical distribution bounded by 80 
    degrees from shore normal. 
    A. Haghighat 'Monte Carlo Methods for Particle Transport' CRC Press 2015'''
    wave_angles = np.load('wave_angles.npy')
    len_wa = 72-1 # Length of wave angle observations minus 1
    i = int(len_wa*rnv1)
    direction = wave_angles[i]+rnv2*(wave_angles[i+1]-wave_angles[i])
    return direction
    
# <codecell>
#%%time
import chaospy as cp
#import pandas as pd
n_samples = 50000*30

joint = cp.J(cp.Uniform(lo=0,up=1),
             cp.Uniform(lo=0,up=1),
             cp.Uniform(lo=0,up=1),
             cp.Uniform(lo=0,up=1),
             cp.Uniform(lo=0,up=1)
            )
clayton = cp.Clayton(joint, theta=2) #why 2?
samples = clayton.sample(size=n_samples)
independent = [cp.Uniform(lo=0,up=1).sample(size=n_samples),
               cp.Uniform(lo=0,up=1).sample(size=n_samples),
               cp.Uniform(lo=0,up=1).sample(size=n_samples),
               cp.Uniform(lo=0,up=1).sample(size=n_samples)
              ]
rnv = {'length':samples[0],
       'hsig':samples[1],
       'a_hsig':samples[2],
       'tps':samples[3],
       'a_tps':samples[4],
       'interim':independent[0], #needs to be reformulated
       'tide':independent[1],
       'direction1':independent[2],
       'direction2':independent[3]
      }

storms = pd.DataFrame()
storms = storms.append({'length':0,
                        'hsig':0,
                        'a_hsig':0,
                        'tps':0,
                        'a_tps':0,
                        'interim':0,
                        'tide':0,
                        'time_end':0,
                        'direction':0},
                        ignore_index=True)

for i in range(n_samples):
    storm = {'length':get_storm_len(rnv['length'][i]),
             'hsig':get_hsig(rnv['hsig'][i]),
             'a_hsig':get_a_hsig(rnv['a_hsig'][i]),                  
             'tps':get_tps(rnv['tps'][i]),                  
             'a_tps':get_a_tps(rnv['a_tps'][i]),
             'interim':get_interim(rnv['interim'][i],storms['time_end'][i]),
             'tide':get_tide(rnv['tide'][i]),
             'direction':get_direction(rnv['direction1'][i],
                                       rnv['direction2'][i])
            }
    storms = storms.append(storm, ignore_index=True)
    storms['time_end'][i+1] = storms['time_end'][i] + \
                              storms['interim'][i] + \
                              storms['length'][i]     
storms = storms.drop(storms.index[[0]])

# <codecell>
#df = pd.DataFrame([storms.hsig,storms.tps,storms.tide,storms.length,storms.interim])
#df = df.T
#pd.tools.plotting.scatter_matrix(df,alpha=0.3)
         
# <codecell>
#rnv = 1-exp(integ)
#storms
   
# <codecell>
#i = 0
#Time = [0]
#max_time = 1000*8766

# <codecell>
#%%timeit
i = 0
Time = [0]
with open(outputfile, 'w') as outfile:
    header = " interim,  length,    tide,    hsig,\
                a_hsig,     tps,   a_tps,   angle\n"
    outfile.write(header)
    for storm in storms.T:
        raw_line = "{:>8},{:>8},{:>8.3f},{:>8.2f},\
                    {:>8.2f},{:>8.2f},{:>8.2f},{:>8.1f}\n"
        line = raw_line.format(storms['interim'][storm],
                               storms['length'][storm],
                               storms['tide'][storm],
                               storms['hsig'][storm],
                               storms['a_hsig'][storm],
                               storms['tps'][storm],
                               storms['a_tps'][storm],
                               storms['direction'][storm])
        outfile.write(line)
        
        Time.append(Time[i]+storms['interim'][storm]+storms['length'][storm])
        i += 1
outfile.close()
pd.to_pickle(storms, 'temp_storms4.npy')
# <codecell>
'''
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
'''