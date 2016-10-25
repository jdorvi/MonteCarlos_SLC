# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 12:25:00 2016

@author: jdorvinen
"""
# <codecell>
import numpy as np
import matplotlib.pyplot as plt
# <codecell>
def brent_solver(func,a,b,tolerance,te,rnv):
    A = func(a,te,rnv)
    B = func(b,te,rnv)
    
    if A*B >= 0:
        return ValueError("Error, root is not bracketed")
    if abs(A) < abs(B):
        (a, b) = (b, a)
    c = a
    s = b
    mflag = 1
    # Solver kernel
    while func(b,te,rnv)!=0 and func(s,te,rnv)!=0 and abs(b-a)>tolerance:
        f_a = func(a,te,rnv)
        f_b = func(b,te,rnv)
        f_c = func(c,te,rnv)
        
        if f_a!=f_c and f_b!=f_c:
            # Inverse quadratic interpolation
            s = (a*f_b*f_c)/((f_a-f_b)*(f_a-f_c)) + \
                (b*f_a*f_c)/((f_b-f_a)*(f_b-f_c)) + \
                (c*f_a*f_b)/((f_c-f_a)*(f_c-f_b))
        else:
            # Secant method
            s = b - f_b*(b-a)/(f_b-f_a)
        
        # Check conditions
        cond1 = (s==min(s,b,((3*a+b)/4)) or s==max(s,b,((3*a+b)/4)))
        if mflag==1:
            cond2 = (abs(s-b)>=abs(b-c)/2)
            cond3 = False
            cond4 = (abs(b-c)<abs(tolerance))
            cond5 = False
        else:
            cond2 = False
            cond3 = (abs(s-b)>=abs(c-d)/2)
            cond4 = False
            cond5 = (abs(c-d)<abs(tolerance))
        # Are any of the conditions true?
        if cond1 or cond2 or cond3 or cond4 or cond5:
            # Bisection method
            s = (a+b)/2
            mflag=1
        else:
            mflag=0
        
        f_s = func(s,te,rnv)
        d = c
        c = b
        if f_a*f_s<0:
            b = s
        else:
            a = s
        if abs(func(a,te,rnv))<abs(func(b,te,rnv)):
            (a, b) = (b, a)
    return (s,rnv)

def func(Gi, te, rnv):
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

# <codecell>
%%time
rnvs = np.arange(0.001,1.00,0.001)
s = []
te = 0.79

for rnv in rnvs:
    s.append(brent_solver(func, 
                          a = 0,
                          b = 1,
                          tolerance = 0.00005,
                          te = te,
                          rnv = rnv))
"""
s_b = []
for rnv in rnvs:
    s_b.append(opt.zeros.brentq(func, 0, 1, 
                                args=(te,rnv),
                                xtol=2e-12,
                                rtol=8.8817841970012523e-16,
                                maxiter=100,
                                full_output=False,
                                disp=True))
"""
# this implementation of the brent solver is ~ 7 times slower than the scipy
# built-in function. Will keep this script for reference but use scipy function
# in final script.
# <codecell>
df1 = pd.DataFrame(s)
df2 = pd.Series(s_b)

df1[0].plot()
df2.plot()