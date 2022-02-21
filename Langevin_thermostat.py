# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 12:49:59 2021

@author: youri
"""

import numpy as np
import matplotlib.pylab as plt


Avogadro = 6.02214086e23
Boltzmann = 1.38064852e-23

def run(**args):
    natoms, temp, mass, radius = args['natoms'], args['temp'], args['mass'], args['radius']
    dt, steps, freq, box = args['dt'], args['steps'], args['freq'],args['box']
    
    positions = np.random.rand(natoms, ndim)
    for i in range(ndim):
        positions[:,i] *= box[i][0] + (box[i][1] - box[i][0])
    
    for x,y,z in positions:
        pass
        

ndim = 3
params = {
    'natoms': 1000,
    'temp': 300,
    'mass': 0.001,
    'radius': 120e-12,
    'relax' : 1e-13,
    'dt' : 1e-15,
    'steps' : 10000,
    'freq' : 100,
    'box' : ((0,1e-8),(0,1e-8),(0,1e-8))
    }
run(**params)