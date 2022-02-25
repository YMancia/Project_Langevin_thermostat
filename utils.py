# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 17:09:27 2022

@author: youri
"""

import numpy as np
import datetime
import matplotlib.pylab as plt

def inputParams():
    """This function creates a parameter dictionary with default values, then uses the values
    from the config file replacing the ones in initial dictionary, if a parameter is not 
    present in the config file, it keeps the default values."""
    try:
        log('---STARTING PARAMETERS LOADING---')

        params = {'natoms': 500, 
                  'temp' : 300, 
                  'mass' : 0.001, 
                  'radius' : 120e-12, 
                  'relax' : 1e-13,
                  'timestep': 1e-15, 
                  'maxsteps' : 3000, 
                  'outputfreq' : 10, 
                  'outputfile': 'tray-Langevin-thermo.dump', 
                  'box' : ((0, 1e-8), (0, 1e-8), (0, 1e-8))}
    
        defaultParams = params.copy()
        with open("config.txt",'r') as file:
            for line in file:
                values = line.split()
                (key, value) = values[0], values[1]
                if key in params.keys() and key != 'box':
                    log(f'INFO\tSetting the parameter {key} = {value}.')
                    params[key] = value
                    defaultParams.pop(key)
            
                elif key == 'box':
                    box = ((float(values[1]),float(values[2])), (float(values[3]),float(values[4])), 
                           (float(values[5]),float(values[6])))
                    params['box'] = box
                    defaultParams.pop('box')
                    log(f'INFO\tSetting the parameter box = {box}.')
            
                else:
                    log(f'WARNING\tNo parameter name found for {key} - skipping it.')
        for key in defaultParams.keys():
            log(f'INFO\tSetting {key} to default value.')
            
    except Exception as err:
        error = str(err)
        log(f'An error occurred while loading the input parameters: {error}')
    finally:
        log('---ENDING PARAMETERS LOADING---')
        return params
    
def ClearOutput(filename):
    """clears the output file"""
    path = './Output/'
    open(path + filename + '.dump', 'w')

def ClearLog():
    """clears the log file"""
    logFile = 'Langevin-simulation-log.txt'
    open(logFile, 'w')
    
def log(string):
    """a simple logging function with date and time"""
    with open('Langevin-simulation-log.txt', 'a') as fp:
        fp.write(str(datetime.datetime.now()) + '\t' + string + '\n')

def writeOutput(filename, natoms, timestep, box, positions, velocities, radius):
    """This function writes the positions and velocities of each particle at each timestep in a
    dump file, this format can be opened with a molecular visualization tool"""
    path = './Output/'
    fp = open(path + filename + '.dump', 'a')
    fp.write('ITEM: TIMESTEP\n')
    fp.write(f'{timestep}\n')
    fp.write('ITEM: NUMBER OF ATOMS\n')
    fp.write(f'{natoms}\n')
    fp.write('ITEM: BOX BOUNDS f f f\n')
    for (a, b) in box:
        fp.write(f'{a} {b}\n')
    fp.write('ITEM: ATOMS radius x y z vx vy vz\n')
    for atom in range(natoms):
        fp.write(f'{radius[atom]} {positions[atom,0]} {positions[atom,1]} {positions[atom,2]} {velocities[atom,0]} {velocities[atom,1]} {velocities[atom,2]}\n')
                
def plot(output):
    """plots the time in picoseconds (1e12) vs the calculated temperature"""
    path = './Output/'
    plt.plot(output[:,0]*1e12, output[:,1])
    plt.xlabel('Time (ps)')
    plt.ylabel('Temp (K)')
    plt.savefig(path + 'Simulation_Temperature.png')
    plt.show()
    
        