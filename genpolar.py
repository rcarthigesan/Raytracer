# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:41:49 2016

@author: raman
"""

import numpy as np
import matplotlib.pyplot as pl

"""
Module containing two generators:
    
    1. rtpairs(R,N): produces a list of polar coordinate pairs, evenly distributed at radii defined in the list R, and with a number of points at each radius defined in the list N.
    
    2. rtuniform(n, rmax, m): produces a list of polar coordinate pairs evenly distributed over a disc of radius rmax in n rings around a central point at (0,0). The ith ring from the centre has m*i evenly distributed points.
    
    Generators include plotting functions
    
"""

def rtpairs(R, N):
    i=0
    for r in R:
        a=(2*(np.pi))/N[i] #angle interval
        i+=1
        for t in range(int(N[i-1])):
            pl.plot(r * np.cos(t*a), r * np.sin(t*a), 'bo')
            yield (r, t*a)
    pl.axis('equal')
    pl.grid()
            
def rtuniform(n, rmax, m):
    R, N=np.zeros(n + 1), np.zeros(n + 1)
    for j in range(int(n+1)):
        R[j]=j*(np.float(rmax)/n)
        N[j]=int(m*j)
    R = np.append(R, 0)
    N = np.append(N, 1)
    i=0
    for r in R:
        a=(2*(np.pi))/N[i] #angle interval
        i+=1
        for t in range(int(N[i-1])):
            #pl.plot(r * np.cos(t*a), r * np.sin(t*a), 'bo')
            yield (r, t*a)
    #pl.grid()
    #pl.axis('equal')