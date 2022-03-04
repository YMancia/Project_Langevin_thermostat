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
import utils
import matplotlib.pylab as plt

def MotionIntegration(positions, velocities, forces, mass, dt):
    """ A simple forward Euler integrator that moves the system in time 
    positions: atomic positions (ndarray, updated)
    velocities: atomic velocity (ndarray, updated)
    """
    positions += velocities * dt
    velocities += forces * dt / mass

def EnforceWallReflection(positions, velocities, box):
    
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
            if (positions[atom,i] >= box[i][1] or positions[atom,i] <= box[i][0]):
                velocities[atom, i] = velocities[atom, i]*(-1)
                
                
def CalculateForces(velocities, relax, mass, temp, timeStep):
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
    forces = np.random.normal(0, sigma, velocities.shape) - (velocities*mass)/relax
    return forces

def ComputeInstTemp(velocities, mass):
    """ Computes the Temperature starting from the velocities of the particles
    @mass: particle mass (ndarray)
    @velocities: particle velocities (ndarray)
    @temp: temperature (float)
    @relax: thermostat constant (float)
    @timestep: simulation timestep (float)
    returns Temperature (float)
    """
    natoms, ndim = velocities.shape
    temp = 0
    for vx, vy, vz in velocities:
        temp += (vx**2 + vy**2 + vz**2)*mass/(boltzmann*natoms*ndim)
    return temp

def CoordInitialitation(box, natoms, ndim):
    """This function initializes the position and velocity arrays with random numbers
    between 0 and 1 with dimension (natoms, ndim). Then scales the positions so that all the
    atoms are contained in the box.
    @box : box size (2, ndim) tuple
    @natoms: number of atoms (int)
    @ndim: number of dimensions (int)
    """
    
    positions = np.random.rand(natoms, ndim)
    velocities = np.random.rand(natoms, ndim)
        
    """Updating the positions in order to have the molecules all inside the box"""
    for i in range(ndim):
        positions[:,i] *= box[i][0] + (box[i][1] - box[i][0])
    return positions, velocities

def RunSimulation(**args):
    """ Takes in input the parameters of a the simulation, calls the initialization function, 
    and calls the force-computing, position-velocities updating and buondary enforcement 
    functions until the max number of steps is reached, returning an array made from timestep
    and temperature calculated at each step.
    @**args: dictionary with the input parameters
    returns output: array with (timestep*step, Temperature) (2,natoms)
    """
    try:
        utils.log('---STARTING SIMULATION---')
        
        output = []
        
        """ Setting up the parameters from the args passed
        @args: dictionary with the input parameters
        """
        
        natoms, temp, mass, radius = args['natoms'], args['temp'], args['mass'], args['radius']
        dt, maxSteps, freq, box, relax = args['timestep'], args['maxsteps'], args['outputfreq'], args['box'] ,args['relax']
        outputFile = args['outputfile']
        
        """setting the mass as the molecular mass"""
        mass = mass/avogadro
        
        """Initializing the velocities and the positions"""
        positions, velocities = CoordInitialitation(box, natoms, ndim)
        

        """starting the iteration"""
        for step in range(maxSteps):
            
            """calculating the forces"""
            forces = CalculateForces(velocities, relax, mass, temp, dt)
            
            """integrating"""
            MotionIntegration(positions, velocities, forces, mass, dt)
            
            """applying reflective boundaries"""
            EnforceWallReflection(positions, velocities, box)
            
            """exporting the output using the output exporting frequency set in the
            config parameters"""
            output.append([dt*step, ComputeInstTemp(velocities, mass)])
            if step%freq == 0:
                utils.writeOutput(outputFile, natoms, dt, box, positions, velocities, np.ones(natoms)*radius)
                utils.log(f'INFO\tExporting output to {outputFile} at step {step}.')
                
                
    except Exception as err:
        error = str(err)
        utils.log(f'An error occurred while the simulation was running: {error}')
    utils.log('---ENDING SIMULATION---')
    return np.array(output)  

"""declaring global constants"""
avogadro = 6.02214086e23
boltzmann = 1.38064852e-23
ndim = 3
"""clearing the log file"""
utils.ClearLog()
"""input the parameters"""
params = utils.inputParams('config.txt')

"""clears the output file"""
utils.ClearOutput(params['outputfile'])

"""runs the simulation"""
output = RunSimulation(**params)

"""plots the output"""
utils.plot(output)

