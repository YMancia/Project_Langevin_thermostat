# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 17:35 2022

@author: youri mancia
"""
# -------------------------------------------------------------------------
#
#   A simple molecular dynamics solver that simulates the motion
#   of non-interacting particles in the canonical ensemble using
#   a Langevin thermostat.
#
# --------------------------------------------------------------------------

import numpy as np
import matplotlib.pylab as plt
import utils

def Integration(positions, velocities, forces, mass, dt):
    """ A simple forward Euler integrator that moves the system in time 
    positions: atomic positions (ndarray, updated)
    velocities: atomic velocity (ndarray, updated)
    """
    positions += velocities * dt
    velocities += forces * dt / mass

def CheckWall(positions, velocities, box):
    
    """ This function enforces reflective boundary conditions.
    All particles that hit a wall have their velocity updated
    in the opposite direction.
    @positions: atomic positions (ndarray)
    @velocities: atomic velocity (ndarray, updated if collisions detected)
    @box: simulation box size (tuple)
    """
    natoms, ndim = positions.shape
    for i in range(ndim):
        for atom in range(natoms):
            if (positions[atom,i] > box[i][1] or positions[atom,i] < box[i][0]):
                velocities[atom, i] = velocities[atom, i]*(-1)
                
                
def CalculateForces(positions, velocities, relax, mass, temp, timeStep):
    """ Computes the Langevin force for all particles
    @mass: particle mass (ndarray)
    @velocities: particle velocities (ndarray)
    @temp: temperature (float)
    @relax: thermostat constant (float)
    @timestep: simulation timestep (float)
    returns forces (ndarray)
    """
    natoms, ndim = velocities.shape
    sigma = np.sqrt(2*boltzmann*temp*mass/(relax*timeStep))
    forces = np.random.normal(0, sigma, positions.shape) - (velocities*mass)/relax
    return forces

def computeInstTemp(velocities, mass):
    """ Computes the Temperature starting from the velocities of the particles
    @mass: particle mass (ndarray)
    @velocities: particle velocities (ndarray)
    @temp: temperature (float)
    @relax: thermostat constant (float)
    @timestep: simulation timestep (float)
    returns forces (ndarray)
    """
    natoms, ndim = velocities.shape
    temp = 0
    for vx, vy, vz in velocities:
        temp += (vx**2 + vy**2 + vz**2)*mass/(boltzmann*natoms*ndim)
    return temp

def run(**args):
    try:
        utils.log('---STARTING SIMULATION---')
        
        nsteps = 0
        output = []
        
        """ Setting up the parameters from the args passed
        @args: dictionary with the input parameters
        """
        
        natoms, temp, mass, radius = int(args['natoms']), float(args['temp']), float(args['mass']), float(args['radius'])
        dt, steps, freq, box, relax = float(args['timestep']), int(args['maxsteps']), int(args['outputfreq']),args['box'],float(args['relax'])
        outputFile = args['outputfile']
        
        """Declaring the positions and the velocities arrays
        with natoms rows and ndim columns
        """
        positions = np.random.rand(natoms, ndim)
        velocities = np.random.rand(natoms, ndim)
        
        """setting the mass as the molecular mass"""
        mass = mass/avogadro
        
        """Updating the positions in order to have the molecules all inside the box"""
        for i in range(ndim):
            positions[:,i] *= box[i][0] + (box[i][1] - box[i][0])

        """starting the iteration"""
        while (nsteps < steps):
            nsteps += 1
            
            """calculating the forces"""
            forces = CalculateForces(positions, velocities, relax, mass, temp, dt)
            
            """integrating"""
            Integration(positions, velocities, forces, mass, dt)
            
            """applying reflective boundaries"""
            CheckWall(positions, velocities, box)
            
            """exporting the output using the output exporting frequency set in the
            config parameters"""
            output.append([dt*nsteps, computeInstTemp(velocities, mass)])
            if nsteps%freq == 0:
                utils.writeOutput(outputFile, natoms, dt, box, positions, velocities, np.ones(natoms)*radius)
                utils.log(f'INFO\tExporting output to {outputFile} at step {nsteps}.')
                
                
    except Exception as err:
        error = str(err)
        utils.log(f'An error occurred while the simulation was running at step {nsteps}: {error}')
    utils.log('---ENDING SIMULATION---')
    return np.array(output)  

"""declaring global constants"""
avogadro = 6.02214086e23
boltzmann = 1.38064852e-23
ndim = 3
"""clearing the log file"""
utils.ClearLog()
"""input the parameters"""
params = utils.inputParams()
utils.ClearOutput(params['outputfile'])
output = run(**params)

"""plots the time in picoseconds (1e12) vs the calculated temperature"""
plt.plot(output[:,0]*1e12, output[:,1])
plt.xlabel('Time (ps)')
plt.ylabel('Temp (K)')
plt.show()
