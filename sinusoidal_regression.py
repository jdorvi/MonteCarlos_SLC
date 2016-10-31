# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 15:02:11 2016

@author: jdorvinen
"""
# <codecell>
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pymc3 as pm
# <codecell>
file_name = 'montauk_combined_data.csv'
df = pd.read_csv(file_name)
df.rename(columns=lambda x: x.strip().rstrip(), inplace=True)
df.end = pd.to_datetime(df.end)
df.set_index('end', inplace=True)
x = df['hsig'].groupby([df.index.month]).agg('count')/13
data = np.matrix((x.index.values,x.values))
data = data.T
dataf = pd.DataFrame(data=data)
#dataf.to_csv('count_month.csv', index=False, header=['month','frequency'])
# <codecell>
# The following code is directly modelled after a blog post by Justin Bozonier
# on using the mcpy3 Python package. 
# Source: http://www.databozo.com/2014/01/17/Exploring_PyMC3.html
trace = None
x_data = dataf[0]
x_data[12] = 13
x_data = (x_data-1)/12
y_data = dataf[1]
y_data[12] = y_data[0]
# <codecell>
def graph(formula, x_range, color='black', alpha=1):
    x = np.array(x_range)
    y = eval(formula)
    plt.plot(x, y, color=color, alpha=alpha)
# <codecell>
with pm.Model() as model:
    alpha = pm.Normal('alpha', mu=0, sd=20)
    beta = pm.Normal('beta', mu=0, sd=20)
    gamma = pm.Normal('gamma', mu=0, sd=20)
    sigma = pm.Uniform('sigma', lower=0, upper=20)
    y_est = alpha + beta*np.sin(2*np.pi*x_data) + gamma*np.cos(2*np.pi*x_data)
    
    likelihood = pm.Normal('y', mu=y_est, sd=sigma, observed=y_data)
    
    start = pm.find_MAP()
    step = pm.NUTS(state=start)
    trace = pm.sample(2000, step, start=start, progressbar=True) #very slow without theano

pm.traceplot(trace)

# <codecell>
plt.scatter(x_data,12*y_data)
for i in np.arange(0,1000):
    #point = trace.point(i)
    formula = '{0} + {1}*np.sin(2*np.pi*{3}) + {2}*np.cos(2*np.pi*{3})'
    #graph(formula.format(point['alpha'],
    #                     point['beta'],
    #                     point['gamma']), 
    #      np.arange(0,13/12,1/12),
    #      color='black',
    #      alpha=0.01)
# Model data fit: alpha=1.498, beta=-0.348, gamma=1.275
x10_data = np.arange(0,1.01,0.01)
y_calc = [eval(formula.format(12*1.498, 12*-0.348, 12*1.275, i)) for i in x10_data]
plt.plot(x10_data, y_calc, color='red')
#1.273,1.106,0.81,7.43)), color='red')
