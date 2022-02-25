# Molecular dynamics: Langevin thermostat
The Langevin thermostat is a molecular simulation tool which maintains the temperature of a system through a modification of the equations of motions:

<img src="https://render.githubusercontent.com/render/math?math=ma = F_{rand} -m\gamma v">

The simulation proceeds through timesteps, at each time step all particles receive a random force and have their velocities lowered using a constant friction. The average magnitude of the random forces and the friction are related in a particular way, which guarantees that the ``fluctuation-dissipation'' theorem is obeyed, thereby guaranteeing NVT statistics. The random contribution of the forces are calculated using the normal distribution with a standard deviation proportional to the temperature:

<img src="https://render.githubusercontent.com/render/math?math=\sigma_F = \sqrt{\frac{2mk_BT}{\gamma \Delta t}}">

where <img src="https://render.githubusercontent.com/render/math?math=\Delta T"> is the timestep of the simulation.

Reflecting boundary conditions are applied to the system: each time a particles finds itself at the box borders its velocity is changed to the opposite direction.

# Simulation steps
1. The simulation starts by acquiring the input parameters, which are the number of atoms, the molecular mass in Kg, the temperature in Kelvin, the radius of the atoms in m, the friction constant ('relax'), the timestep, the max number of steps the 3D box size and the output values such as the frequency of output update and the output filename.
2. The positions and the initial velocities of the atoms are created using an uniform distribution, then the positions are re-scaled to the box size.
3. At each step, an array of forces is calculated and the positions and the velocities are update through a forward Euler integrator. Then the updated positions are checked and if any wall has been hit, the direction of the velocity for which the collision was made is multiplied by -1. 
4. At the end of the simulation a .dump file is created in order to store the trajectories of the particles. This file can be opened by using a simple visualization tool (like OVITO) and the whole timeline can be visualized.
5. The whole simulation is tracked by a log which provides any error, warning or info to the user.

#Structure of the project
1. The config.txt file has to be written in order to setup the initial parameters. If any of the parameters is not provided in the file the software will use the default values, informing the user through the simulation log ('Langevin-simulation-log.txt').
2. The module 'utils.py' contains the input and output, the logging and the plotting functions.
3. The main simulation is contained in the 'Langevin_thermostat.py' module


#Installation
The project requires: numpy, datetime and matplotlib.
Install the repository and the requirements by typing:
```
git clone https://github.com/Mitenus/Project_Langevin_thermostat
pip install numpy
pip install datetime
pip install matlplotlib
```
In order to visualize the simulation trajectory you will need to install a visualization tool, OVITO is suggested: https://www.ovito.org/
Open the tool and import the data as LAMMPS text dump file ('auto-detect format file' would cause a crash). Set the visualization speed at 60 fps for a clear view.

