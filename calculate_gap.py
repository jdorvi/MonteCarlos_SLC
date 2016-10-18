# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 17:57:40 2016

@author: jdorvinen
"""
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# <codecell>
# Model data fit: alpha=1.498, beta=-0.348, gamma=1.275
# Callaghan et al. used: alpha=21.46, beta=1.08, gamma=1.07
a = 12*1.498
b = 12*-0.348
c = 12*1.275
w = 2*np.pi
rnv = np.random.random()

# Takes a random variable and can be used to find a value for Gi
formulaG = '1 - np.exp(-(a*w*Gi \
            + b*(np.cos(w*te) - np.cos(w*(te + Gi))) \
            - c*(np.sin(w*te) - np.sin(w*(te + Gi))))/w)'

# Initial estimate of Gi. Obtained from the second order Taylor series 
# expansion about Gi=0 of "formulaG"
formulaGi_0 = 'rnv / (a + b*np.sin(w*te[i-1]) + c*np.cos(w*te[i-1]))'

def func(te,Gi):
    z = eval(formulaG)
    return z

te = np.arange(0,1.01,0.01)
Gi = np.arange(0,1.01,0.01)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(te, Gi)
zs = np.array([func(te,Gi) for te,Gi in zip(np.ravel(X), np.ravel(Y))])
Z = zs.reshape(X.shape)

zs2 = np.array([0 for te,Gi in zip(np.ravel(X), np.ravel(Y))])
Z2 = zs2.reshape(X.shape)

ax.plot_surface(X,Y,Z, color='blue', alpha=0.5)
ax.plot_surface(X,Y,Z2, color='yellow', alpha=0.8)
ax.set_xlabel('TimeEnd')
ax.set_ylabel('Gi')
ax.set_zlabel('RNV')

plt.show()
        
