# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 13:49:00 2022

@author: youri
"""
import Langevin_thermostat
import utils
import hypothesis
import numpy as np
from hypothesis import strategies as st
from hypothesis import settings
from hypothesis import given

def test_EnforceWallReflection1():
    """This function tests the EnforceWallReflection function by creating a position array
    of 4 atoms at the borders of the box and checking that the velocity after the 
    EnforceWallReflection are reverted in the direction of the hitted wall
    """
    natoms = 4
    ndim = 3
    velocities = np.random.rand(natoms, ndim)
    afterCheckVelocities = velocities.copy()
    box = ((-1e-8, 1e-8), (-1e-8, 1e-8), (-1e-8, 1e-8))
    positions = np.array([[0, 1e-8, 0],
                          [1e-8, 0, 0],
                          [0, 0, -1e-8],
                          [0, 0, 0]])
    Langevin_thermostat.EnforceWallReflection(positions, afterCheckVelocities, box)
    
    assert velocities[0][0] == afterCheckVelocities[0][0]
    assert velocities[0][1] == -1*afterCheckVelocities[0][1]
    assert velocities[0][2] == afterCheckVelocities[0][2]
    assert velocities[1][0] == -1*afterCheckVelocities[1][0]
    assert velocities[1][1] == afterCheckVelocities[1][1]
    assert velocities[1][2] == afterCheckVelocities[1][2]
    assert velocities[2][0] == afterCheckVelocities[2][0]
    assert velocities[2][1] == afterCheckVelocities[2][1]
    assert velocities[2][2] == -1*afterCheckVelocities[2][2]
    assert velocities[3][0] == afterCheckVelocities[3][0]
    assert velocities[3][1] == afterCheckVelocities[3][1]
    assert velocities[3][2] == afterCheckVelocities[3][2]

@given(natoms=st.integers(1,10000)) 
def test_EnforceWallReflection(natoms):
    """This function tests the EnforceWallReflection function by creating a position array of natoms
    at the borders of the box and checking that the velocity after the CheckWall are all
    reverted
    @natoms: number of atoms (int)
    """
    ndim = 3
    box = ((0, 1e-8), (0, 1e-8), (0, 1e-8))
    positionsAtBound = np.ones((natoms, ndim))
    positionsAtCenter = np.ones((natoms, ndim))
    velocities = np.random.rand(natoms, ndim)
    afterCheckVelocitiesBound = velocities.copy()
    afterCheckVelocitiesCenter = velocities.copy()
    for i in range(ndim):
        positionsAtBound[:,i] *= box[i][0] + (box[i][1] - box[i][0])
        positionsAtCenter[:,i] *= box[i][0] + (box[i][1] - box[i][0])/2
    Langevin_thermostat.EnforceWallReflection(positionsAtBound, afterCheckVelocitiesBound, box)
    Langevin_thermostat.EnforceWallReflection(positionsAtCenter, afterCheckVelocitiesCenter, box)
    assert np.array_equal(velocities, -1*afterCheckVelocitiesBound)
    assert np.array_equal(velocities, afterCheckVelocitiesCenter)
    
def test_MotionIntegration1():
    """This function test the MotionIntegration function by checking that if the 
    timestep is 0 the resulting positions and velocities remain unchanged
    """
    mass = 0.001
    ndim = 3
    dt = 0.
    box = ((0, 1e-8), (0, 1e-8), (0, 1e-8))
    positions = np.random.rand(4, ndim)
    velocities = np.random.rand(4, ndim)
    forces = np.random.rand(4, ndim)
    afterIntegVelocities = velocities.copy()
    for i in range(ndim):
        positions[:,i] *= box[i][0] + (box[i][1] - box[i][0])
    afterIntegPositions = positions.copy()
    Langevin_thermostat.MotionIntegration(afterIntegPositions, afterIntegVelocities, forces, mass, dt)
    assert np.array_equal(positions, afterIntegPositions)
    assert np.array_equal(velocities, afterIntegVelocities)
    
def test_MotionIntegration2():
    """This function test the MotionIntegration function by checking that if the 
    function is called two times with opposite timesteps the positions and 
    the velocities remain unchanged
    """
    mass = 0.001
    dt = 1e-15
    ndim = 3
    box = ((0, 1e-8), (0, 1e-8), (0, 1e-8))
    positions = np.random.rand(4, ndim)
    afterIntegPositions = positions.copy()
    velocities = np.random.rand(4, ndim)
    forces = np.random.rand(4, ndim)
    afterIntegVelocities = velocities.copy()
    Langevin_thermostat.MotionIntegration(afterIntegPositions, afterIntegVelocities, forces, mass, dt)
    Langevin_thermostat.MotionIntegration(afterIntegPositions, afterIntegVelocities, forces, mass, -1*dt)
    assert np.array_equal(positions, afterIntegPositions)
    assert np.array_equal(velocities, afterIntegVelocities)

@given(natoms = st.integers(1,10000))
def test_ComputeInstTemp1(natoms):
    """This function test the ComputeInstTemp function by checking that for any atom number
    if the velocities are all zeros, the temperature is zero"""
    ndim = 3
    mass = 0.001
    zeroTemp = Langevin_thermostat.ComputeInstTemp(np.zeros((natoms, ndim)), mass)
    assert zeroTemp == 0

@given(natoms = st.integers(1,10000))
def test_ComputeInstTemp2(natoms):
    """This function test the ComputeInstTemp function by checking that for opposite
    velocities, the temperature remains unchanged"""
    ndim = 3
    velocities = np.random.rand(natoms, ndim)
    mass = 0.001
    Temp = Langevin_thermostat.ComputeInstTemp(velocities, mass)
    reverseVelTemp = Langevin_thermostat.ComputeInstTemp(-1*velocities, mass)
    assert Temp == reverseVelTemp

def test_CoordInitialization():
    """This function test the CoordInitialization function by checking that every 
    atom is placed inside the box"""
    ndim = 3
    natoms = 10
    box = ((0, 1e-8), (0, 1e-8), (0, 1e-8))
    positions, velocities = Langevin_thermostat.CoordInitialitation(box, natoms, ndim)
    for dim in range(ndim):
        for atom in range(natoms):
            assert positions[atom,dim] <= box[dim][1] and positions[atom,dim] >= box[dim][0]
        
def test_inputParams():
    """This function test the inputParams function by checking that:
    for an empty input file, the parameters remain at the default value
    for a custom input file the parameters get correctly modified
    for a differently specified box the box parameter is correclty generated
    """
    defaultParams = utils.inputParams('Config_Test/empty.txt')
    assert defaultParams == {'natoms': 500, 
                             'temp' : 300, 
                             'mass' : 0.001, 
                             'radius' : 120e-12, 
                             'relax' : 1e-13,
                             'timestep': 1e-15, 
                             'maxsteps' : 3000, 
                             'outputfreq' : 10, 
                             'outputfile': 'tray-Langevin-thermo', 
                             'box' : ((0, 1e-8), (0, 1e-8), (0, 1e-8))}
    
    params1 = utils.inputParams('Config_Test/params1.txt')
    assert params1 == {'natoms': 5000, 
                       'temp' : 500, 
                       'mass' : 15, 
                       'radius' : 432, 
                       'relax' : 667,
                       'timestep': 987, 
                       'maxsteps' : 10000, 
                       'outputfreq' : 1000, 
                       'outputfile': 'tray-Langevin-thermo', 
                       'box' : ((0, 1e-8), (0, 1e-8), (0, 1e-8))}
    
    diffBoxParams = utils.inputParams('Config_Test/diffBoxParams.txt')
    assert diffBoxParams == {'natoms': 500, 
                             'temp' : 300, 
                             'mass' : 0.001, 
                             'radius' : 120e-12, 
                             'relax' : 1e-13,
                             'timestep': 1e-15, 
                             'maxsteps' : 3000, 
                             'outputfreq' : 10, 
                             'outputfile': 'tray-Langevin-thermo', 
                             'box' : ((10, 20), (30, 40), (50, 60))}
