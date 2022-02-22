# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 17:35 2022

@author: youri
"""

import numpy as np
import matplotlib.pylab as plt


def Verlet(positions, velocities, forces, mass, dt):
    positions[:,:] = positions[:,:] + velocities[:,:]*dt
    velocities[:,:] = velocities[:,:] + dt*forces[:,:]*(1/mass)

def CheckWall(positions, velocities, box, ndim, natoms):
    for i in range(ndim):
        for atom in range(natoms):
            if (positions[atom,i] > box[i][1] or positions[atom,i] < box[i][0]):
                velocities[atom, i] = velocities[atom, i]*(-1)
                
                
def CalculateForces(positions, velocities, relax, mass, temp, dt):
    sigma = np.sqrt(2*temp*boltzmann*mass/(relax*dt))
    forces = np.random.normal(0, sigma, positions.shape)
    forces[:,:] = forces[:,:] - mass*velocities/relax
    return forces

def computeInstTemp(velocities, mass, natoms):
    temp = 0
    for vx, vy, vz in velocities:
        temp = temp + mass*(vx**2 + vy**2 + vz**2)
    return temp/natoms

def run(**args):
    natoms, temp, mass, radius = args['natoms'], args['temp'], args['mass'], args['radius']
    dt, steps, freq, box, relax = args['dt'], args['steps'], args['freq'],args['box'],args['relax']
    
    output = []
    positions = np.random.rand(natoms, ndim)
    velocities = np.random.rand(natoms, ndim)
    for i in range(ndim):
        positions[:,i] *= box[i][0] + (box[i][1] - box[i][0])
    
    nsteps = 0
    while (nsteps < steps):
        nsteps += 1
        forces = CalculateForces(positions, velocities, relax, mass, temp, dt)
        Verlet(positions, velocities, forces, mass, dt)
        CheckWall(positions, velocities, box, ndim, natoms)
        output.append([dt*nsteps, computeInstTemp(velocities, mass, natoms)])
    return np.array(output)        
    
avogadro = 6.02214086e23
boltzmann = 1.38064852e-23
ndim = 3
params = {
    'natoms': 1000,
    'temp': 300,
    'mass': 0.001,
    'radius': 120e-12,
    'relax' : 1e-13,
    'dt' : 1e-15,
    'steps' : 3000,
    'freq' : 100,
    'box' : ((0,1e-8),(0,1e-8),(0,1e-8))
    }
output = run(**params)
plt.plot(output[:,0]*1e12, output[:,1])
plt.xlabel('Time (ps)')
plt.ylabel('Temp (K)')
plt.show()
