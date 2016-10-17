# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 17:57:40 2016

@author: jdorvinen
"""
import numpy as np

# <codecell>
# Model data fit: alpha=1.498, beta=-0.348, gamma=1.275
a = 1.498
b = -0.348
g = 1.275
w = 2*np.pi
rnv = np.random.random()

# Takes a random variable and can be used to find a value for Gi
formulaG = '1 - np.exp((a*w*Gi \
            + b*(np.cos(w*te[i-1]) - np.cos(w*(te[i-1] + Gi))) \
            - c*(np.sin(w*te[i-1]) - np.sin(w*(te[i-1] + Gi))))/w)'

# Initial estimate of Gi. Obtained from the second order Taylor series 
# expansion about Gi=0 of "formulaG"
formulaGi_0 = 'rnv / (a + b*np.sin(w*te[i-1]) + c*np.cos(w*te[i-1]))'

