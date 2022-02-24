# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 17:35 2022

@author: youri
"""

import numpy as np
import matplotlib.pylab as plt
import utils

def inputParams():

    params = {'natoms': 500, 'temp' : 300, 'mass' : 0.001, 'radius' : 120e-12, 'relax' : 1e-13,
              'timestep': 1e-15, 'maxsteps' : 3000, 'outputfreq' : 10, 'outputfile': 'tray-Langevin-thermo.dump'}
    with open("config.txt") as file:
        for line in file:
            values = line.split()
            (key, value) = values[0], values[1]
            if key in params.keys():
                params[key] = value
            if key == 'box':
                box = ((float(values[1]),float(values[2])), (float(values[3]),float(values[4])), 
                (float(values[5]),float(values[6])))
                params['box'] = box
        return params

def Integration(positions, velocities, forces, mass, dt):
    positions += velocities * dt
    velocities += forces * dt / mass

def CheckWall(positions, velocities, box, ndim, natoms):
    for i in range(ndim):
        for atom in range(natoms):
            if (positions[atom,i] > box[i][1] or positions[atom,i] < box[i][0]):
                velocities[atom, i] = velocities[atom, i]*(-1)
                
                
def CalculateForces(positions, velocities, relax, mass, temp, dt, natoms, ndim):
    sigma = np.sqrt(2*boltzmann*temp*mass/(relax*dt))
    forces = np.random.normal(0, sigma, positions.shape) - (velocities*mass)/relax
    return forces

def computeInstTemp(velocities, mass, natoms, ndim):
    temp = 0
    for vx, vy, vz in velocities:
        temp += (vx**2 + vy**2 + vz**2)*mass/(boltzmann*natoms*ndim)
    return temp

def run(**args):
    natoms, temp, mass, radius = int(args['natoms']), float(args['temp']), float(args['mass']), float(args['radius'])
    dt, steps, freq, box, relax = float(args['timestep']), int(args['maxsteps']), args['outputfreq'],args['box'],float(args['relax'])
    outputFile = args['outputfile']
    output = []
    positions = np.random.rand(natoms, ndim)
    velocities = np.random.rand(natoms, ndim)
    mass = mass/avogadro
    for i in range(ndim):
        positions[:,i] *= box[i][0] + (box[i][1] - box[i][0])

    nsteps = 0
    open(outputFile, 'w')
    while (nsteps < steps):
        nsteps += 1
        forces = CalculateForces(positions, velocities, relax, mass, temp, dt, natoms, ndim)
        Integration(positions, velocities, forces, mass, dt)
        CheckWall(positions, velocities, box, ndim, natoms)
        output.append([dt*nsteps, computeInstTemp(velocities, mass, natoms, ndim)])
        if nsteps%freq == 0:
            utils.writeOutput(outputFile, natoms, dt, box, positions, velocities, np.ones(natoms)*radius)
    return np.array(output)  
          
avogadro = 6.02214086e23
boltzmann = 1.38064852e-23
ndim = 3
params = inputParams()
output = run(**params)
plt.plot(output[:,0]*1e12, output[:,1])
plt.xlabel('Time (ps)')
plt.ylabel('Temp (K)')
plt.show()
