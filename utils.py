# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 17:09:27 2022

@author: youri
"""

import numpy as np

def writeOutput(filename, natoms, timestep, box, positions, velocities, radius):
    
    fp = open(filename, 'a')
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
                
                
        