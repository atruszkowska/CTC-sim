import subprocess, os
import numpy as np
import scipy
from scipy.stats import gaussian_kde
from scipy.optimize import differential_evolution
import shutil
from itertools import product
from collections import defaultdict
import time
import math

import sys
py_src_path = '../../scripts/'
sys.path.insert(0, py_src_path)

import utils as ut
import simulation_utils as smut

# Script for running the optimization simulations
# Goal: find bond properties that best match the experimental data 

def simulation_driver(new_tag, new_value):

    # Settings 
    exe_name         = 'lns'
    compile_script   = 'compilation.py'
    param_file       = 'parameters.txt'
    log_name         = 'sim_progress.log'
    Np = 300 
    Lx = 0.12
    Ly = 1e-3 
    max_area = Lx*Ly
    max_len = 4.0400e-05
    alpha = None
    # First time step
    t0 = '0' 
    # Last time step
    tf = '36000'
    # Time step for starting the averaging
    t1 = '9000'

    # List of WSS values
    wss_values       = [0, 2, 5, 20]

    # Lift coefficient
    lift_coef = {'0' : 0, '2': 0, '5': 0, '20': 0}

    # Error for each shear stress
    rel_error = []
    all_pdf_errors = []
    
    # Experimental data
    exp_data = {}
    with open('experimental_area.txt', 'r') as value_file:
        for key, line in zip(wss_values, value_file):
            value = float(line.strip())
            exp_data[key] = value/100

    # Run the simulations for each shear stress level 
    for wss in wss_values:

        # Shear rate dependent values
        eta_i = 0.27 
        cl_cluster = lift_coef[str(wss)]

        # Substitute the parameter that is changing
        params = {
          'Bond normal stiffness': 2e-5,
          'Bond shear stiffness':  1e-17,
          'Bond normal strength':  5e-10,
          'Bond shear strength':   8e-24,
          'Wall shear stress': wss,
          'Amplitude of noise': eta_i,
          'Single cell lift coefficient': -0.2,
          'Cell cluster lift coefficient': -0.1,
          'Number of particles': 300
        }
        params[new_tag] = new_value 
        ut.sub_parameters(params, param_file)

        # Remove this !!!!!!!!!!!!!!!!!!!!!!
        params['Cell cluster lift coefficient'] = params['Single cell lift coefficient']/2


        # Run the simulation
        subprocess.call(['./' + exe_name], shell = True)

        # Compute the error for this shear stress
        pos_file = 'results' + str(wss) + '.txt'    
        bond_file = 'bonds' + str(wss) + '.txt'    
        
        # Check if both files exist, skip if not
        if not (os.path.exists(pos_file) and os.path.exists(bond_file)):
            continue
    
        # Load data 
        try:
            df_pos = smut.read_pos(pos_file)
            df_bond = smut.read_bonds(bond_file)
            # Detect clusters and compute average areas at each collected time step
            uf_areas, uf_clusters, uf_counts = smut.cluster_data_union_find(df_bond, df_pos, Np, max_area, max_len, alpha)
        except:
            # Write current values to file
            with open(log_name, 'a') as fout:
                fout.write(" ".join(map(str, x)) + f" {wss} \n")
            return None

        # Group data by realization number
        realization_data = defaultdict(dict)
        for key, value in uf_areas.items():
            time_step, ri = key.split('-')
            realization_data[ri][time_step] = value

        # Compute the ratio for each realization
        ratios = []
        
        for ri, times in realization_data.items():
            if t0 in times:
                values = [
                    value for t, value in times.items()
                    if (int(t) >= int(t1)) and not (isinstance(value, float) and math.isnan(value))]        
                if values and not math.isclose(times[t0], 0.0, rel_tol=1e-9, abs_tol=1e-12):
                    average_tf = sum(values) / len(values)
                    ratio = average_tf / times[t0]
                    ratios.append(ratio)
        
        average_ratio = sum(ratios) / len(ratios) if ratios else float('nan')

        # Compare to experimental data for this shear level
        reference = exp_data[wss]
        rel_error.append(abs((average_ratio - reference) / reference)*100 if reference != 0 else float('inf'))

        # Write current values to file
        with open(log_name, 'a') as fout:
            fout.write(new_tag + str(new_value) + f" {wss} {average_ratio} {rel_error[-1]}\n")

        # Rename the files with results (atomic)
        tag_file = new_tag.replace(' ', '_')
        
        for basename in ['results', 'bonds', 'bond_failures']:
            old_name = basename + str(wss) + '.txt'
            new_name = basename + str(wss) + '_' + tag_file + '_' + str(new_value) + '.txt'
        
            os.rename(old_name, new_name)

    # Write the current average errors to file (cluster area error, positioning error, average of the two)
    average_error = sum(rel_error)/len(rel_error)
    with open(log_name, 'a') as fout:
        fout.write(f"Average error: {average_error}\n")

    return average_error 

# Compile the code 
py_version       = 'python3'
subprocess.call([py_version + ' compilation.py'], shell=True)

# Parameters (changing one at a time)
params = {
#    'Bond normal stiffness': [1.4e-5, 2.6e-5],
#    'Bond shear stiffness':  [0.7e-17, 1.3e-17],
#    'Bond normal strength':  [3.5e-10, 6.5e-10],
#    'Bond shear strength':   [5.6e-24, 10.4e-24],
#    'Single cell lift coefficient': [-0.14, -0.26],
    'Single cell lift coefficient': [0, 0.2],
#    'Cell cluster lift coefficient': [-0.143, -0.077],
#     'Number of particles': [200] #[200, 400] 
}

# Run all, one at a time
for key, values in params.items():
    for value in values:
        error = simulation_driver(key, value)
        print(error)    

