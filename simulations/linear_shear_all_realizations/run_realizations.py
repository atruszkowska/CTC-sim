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

def simulation_driver(x):

    # Settings 
    exe_name         = 'arel'
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

    # New parameter set
    kn, ks, sigma, tau, cl_cluster, cl_factor = x
    # Lift coefficient
    lift_coef = {'0' : 0, '2': cl_cluster, '5': cl_cluster, '20': cl_cluster}

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

        # Substitute the new guesses
        params = {
          'Bond normal stiffness': kn,
          'Bond shear stiffness':  ks,
          'Bond normal strength':  sigma,
          'Bond shear strength':   tau,
          'Wall shear stress': wss,
          'Amplitude of noise': eta_i,
          'Single cell lift coefficient': cl_cluster,
          'Cell cluster lift coefficient': cl_cluster/cl_factor
        }
        ut.sub_parameters(params, param_file)

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
                    print(values/times[t0])
                    average_tf = sum(values) / len(values)
                    ratio = average_tf / times[t0]
                    ratios.append(ratio)
        
        average_ratio = sum(ratios) / len(ratios) if ratios else float('nan')

        # Compare to experimental data for this shear level
        reference = exp_data[wss]
        rel_error.append(abs((average_ratio - reference) / reference)*100 if reference != 0 else float('inf'))

        # Write current values to file
        with open(log_name, 'a') as fout:
            fout.write(" ".join(map(str, x)) + f" {wss} {average_ratio} {rel_error[-1]} \n")

    # Write the current average errors to file (cluster area error, positioning error, average of the two)
    average_error = sum(rel_error)/len(rel_error)
    with open(log_name, 'a') as fout:
        fout.write(f"Average error: {average_error}\n")

    return average_error 

# Compile the code 
py_version       = 'python3'
subprocess.call([py_version + ' compilation.py'], shell=True)

# Parameters 
kn_list          = [2e-5] 
ks_list          = [1e-17]
sigma_max_list   = [5e-10] 
tau_max_list     = [8e-24]
cl_cluster = [-0.2]
cl_factor = [2]

configurations = list(product(kn_list, ks_list, sigma_max_list, tau_max_list, cl_cluster, cl_factor))

for kn, ks, sigma, tau, cl_cluster, cl_factor in configurations:
    error = simulation_driver([kn, ks, sigma, tau, cl_cluster, cl_factor])
    print(error)    

