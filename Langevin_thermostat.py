# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 17:35 2022

@author: youri
"""

import numpy as np
import matplotlib.pylab as plt
import utils

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
    try:
        utils.log('---STARTING SIMULATION---')
        nsteps = 0
        natoms, temp, mass, radius = int(args['natoms']), float(args['temp']), float(args['mass']), float(args['radius'])
        dt, steps, freq, box, relax = float(args['timestep']), int(args['maxsteps']), int(args['outputfreq']),args['box'],float(args['relax'])
        outputFile = args['outputfile']
        output = []
        positions = np.random.rand(natoms, ndim)
        velocities = np.random.rand(natoms, ndim)
        mass = mass/avogadro
        for i in range(ndim):
            positions[:,i] *= box[i][0] + (box[i][1] - box[i][0])

        
        while (nsteps < steps):
            nsteps += 1
            forces = CalculateForces(positions, velocities, relax, mass, temp, dt, natoms, ndim)
            Integration(positions, velocities, forces, mass, dt)
            CheckWall(positions, velocities, box, ndim, natoms)
            output.append([dt*nsteps, computeInstTemp(velocities, mass, natoms, ndim)])
            if nsteps%freq == 0:
                utils.writeOutput(outputFile, natoms, dt, box, positions, velocities, np.ones(natoms)*radius)
                utils.log(f'INFO\tExporting output to {outputFile} at step {nsteps}.')
    except Exception as err:
        error = str(err)
        utils.log(f'An error occurred while the simulation was running at step {nsteps}: {error}')
    utils.log('---ENDING SIMULATION---')
    return np.array(output)  
          
avogadro = 6.02214086e23
boltzmann = 1.38064852e-23
ndim = 3
utils.ClearLog()
params = utils.inputParams()
utils.ClearOutput(params['outputfile'])
output = run(**params)
plt.plot(output[:,0]*1e12, output[:,1])
plt.xlabel('Time (ps)')
plt.ylabel('Temp (K)')
plt.show()
