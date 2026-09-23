import os
import sys
import subprocess
import time
from itertools import product
import numpy as np

sys.path.insert(0, '../../../scripts/')
import utils as ut

### Configuration 
py_version       = 'python3'
exe_name         = 'p_exe'
compile_script   = 'compilation.py'
param_file       = 'parameters.txt'
outfile = "study_results_parametric.txt"

### Inputs
filename = "bonds.txt"
pos_file = "results.txt"        
# Number of cells
N = 5   
last_timestep = 1200
duration = 120
max_realization = 0
# Scaling factor for y-coordinate
scaling_factor = 8.51e-6         

# Parameters to be tested 
sweep = {}
sweep['Wall bond normal stiffness'] = [1e-9] 
sweep['Wall bond shear stiffness'] = [1e-10] 
sweep['Wall bond normal strength'] = [1e-14] 
sweep['Wall bond shear strength'] = [1e-16] 
sweep['Wall bond max gap'] = [6*scaling_factor] 

final_parameters = dict.fromkeys(sweep.keys(), 0.0)

# Generate all combinations from sweep values
configurations = list(product(*sweep.values()))

# Simulations
for kn, ks, sigma, tau, gap in configurations:

    ### Substitute
    final_parameters['Wall bond normal stiffness'] = kn
    final_parameters['Wall bond shear stiffness'] = ks
    final_parameters['Wall bond normal strength'] = sigma
    final_parameters['Wall bond shear strength'] = tau
    final_parameters['Wall bond max gap'] = gap
   
	### Run the simulation
    # Update parameters.txt
    ut.sub_parameters(final_parameters, param_file)
    # Compile 
    ut.subprocess.run([py_version + ' ' + compile_script], shell=True)
    # Run        
    ut.subprocess.call(['./' + exe_name], shell=True)

    ### Compute the metrics
    # Check if all cells are still in the cluster
    # 0 - not in, 1 - in, average over all realizations
    
    # Results: one flag per realization
    results = np.zeros(max_realization + 1, dtype=int)

    # Track cell appearances per realization
    seen = [set() for _ in range(max_realization + 1)]

    with open(filename) as f:
        for line in f:
            parts = line.strip().split()
            cellpair, timestep, realization = parts
            timestep = int(timestep)
            realization = int(realization)

            if timestep != last_timestep:
                continue

            c1, c2 = map(int, cellpair.split('-'))
            seen[realization].update([c1, c2])

    # Now decide per realization
    for r in range(max_realization + 1):
        results[r] = 1 if len(seen[r]) == N else 0

    # Compute average
    average_fraction = results.mean()

    # Compute how far the outermost cell is from the wall
    # y-position of outermost cell per realization
    outermost_y = np.zeros(max_realization + 1)
    with open(pos_file) as f:
        for line in f:
            parts = line.strip().split()
            x, y, timestep, realization = parts
            x, y = float(x), float(y)
            timestep, realization = int(timestep), int(realization)
            
            if timestep != last_timestep:
                continue
            
            # Update if this y is larger than current
            if y > outermost_y[realization]:
                outermost_y[realization] = y    

    # Compute average y, scaled
    average_outermost_y = outermost_y.mean() / scaling_factor

    # Compute cluster velocity
    # Dictionary to store x positions per cell and realization
    # { realization : { cell_index : [x0, x_last] } }
    cell_positions = {r: {} for r in range(max_realization + 1)}
    # Read file and record x positions at t=0 and t=last_timestep
    with open(pos_file) as f:
        for line in f:
            parts = line.strip().split()
            x, y, timestep, realization = parts
            x = float(x)
            timestep = int(timestep)
            realization = int(realization)

            # Assign cell index per line (assume rows correspond to cells 0..N-1 in order)
            if timestep == 0:
                for i in range(0, N):
                    parts = line.strip().split()
                    x, y, timestep, realization = parts
                    x = float(x)
                    timestep = int(timestep)
                    realization = int(realization)
                    cell_positions[realization][i] = [x, None]
                    if not (i == N-1):  
                        line = next(f)

            if timestep == last_timestep:
                for i in range(0, N):
                    parts = line.strip().split()
                    x, y, timestep, realization = parts
                    x = float(x)
                    timestep = int(timestep)
                    realization = int(realization)
                    cell_positions[realization][i][1] = x
                    if not (i == N-1):  
                        line = next(f)

    # Compute average velocity per realization
    velocities = np.zeros(max_realization + 1)

    for r in range(max_realization + 1):
        v_cells = []
        for cell, pos in cell_positions[r].items():
            x0, x_last = pos
            v = (x_last - x0) / duration 
            v_cells.append(v)
        velocities[r] = np.mean(v_cells)

    # Average over all realizations
    average_velocity_all_realizations = velocities.mean()

    ### Append to file - parameter set and metrics
    with open(outfile, "a") as f:
        params_str = (
            f"{final_parameters['Wall bond normal stiffness']} "
            f"{final_parameters['Wall bond shear stiffness']} "
            f"{final_parameters['Wall bond normal strength']} "
            f"{final_parameters['Wall bond shear strength']} "
            f"{final_parameters['Wall bond max gap']} "
        )
        f.write(f"Parameters: {params_str}\n")
        f.write(f"Cells in the cluster: {average_fraction}\n")
        f.write(f"Average outermost cell y position: {average_outermost_y}\n")
        f.write(f"Average velocity over all realizations: {average_velocity_all_realizations}\n")


