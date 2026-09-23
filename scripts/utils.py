# # # # # # # # # # # # # # # # # #
# Utility functions for Python code
# # # # # # # # # # # # # # # # # #

import sys, os, subprocess
import time
import shutil
from colors import *

def msg(text, clr):
    ''' Print text with color defined as clr and style font '''
    print(clr + text + RESET)

def test_pass(res, name):
    ''' Print results of running a test '''
    if res:
        msg('  - ' +  name + ' test passed', GREEN)
    else:
        msg('  ' + name + ' test failed', MAGENTA)

def float_equality(num1, num2, tol):
    ''' Compare two floating point numbers for approximate equality '''
    
    max_num = max([1.0, abs(num1), abs(num2)])
    return abs(num1 - num2) <= tol*max_num
        

def sub_parameters(params_to_sub, fname):
    ''' Substitute specific parameter values 
       
        Assumes format // parameter name in one line and 
            the value in the next line. '''
    
    # Copy the file, substituting where requested
    with open(fname, 'r') as fin, open('temp.txt', 'w') as fout:
        for line in fin:
            if '//' in line:
                # This line contains parameter name
                # Write this line to file (parameter name)
                fout.write(line)
                # First element is //, then its the name
                param_name = (' ').join(line.strip().split()[1:])

                if param_name in params_to_sub:
                    # This parameter needs to be changed
                    # Value should be on the next line  
                    fout.write(str(params_to_sub[param_name]) +'\n')
                    next(fin)
            else:
                fout.write(line)    

    # Move temp to fname        
    shutil.move('temp.txt', fname)

def run_simulations_screens(kn, ks, sigma, tau, wss, param_file, main_name, py_version, compile_script, exe_name):
  
    # Root directory name
    root_dir = os.path.abspath(os.getcwd())
 
    # Directory with these properties
    dir_name = main_name + ('_').join(list(map(str, [kn, ks, sigma, tau, wss])))
    print('Processing ' + dir_name) 

    # Make the directory
    subprocess.run('mkdir ' + dir_name, shell=True, check=True)
    
    # Copy all the files from templates
    subprocess.run('cp -r templates/*.* ' + dir_name + '/', shell=True, check=True)
    
    # Go to directory
    os.chdir(dir_name)

    # Update parameters.txt
    params = {
      'Bond normal stiffness': kn,
      'Bond shear stiffness':  ks,
      'Bond normal strength':  sigma,
      'Bond shear strength':   tau,
      'Wall shear stress': wss
    }
    sub_parameters(params, param_file)
            
    # Compile 
    subprocess.run([py_version + ' ' + compile_script], shell=True)

    # Run in a screen       
    subprocess.call(['screen -d -m ./' + exe_name], shell=True)
                
    # Go back to root
    os.chdir(root_dir) 

def run_simulations_screens_dict(properties, param_file, main_name, py_version, compile_script, exe_name):
  
    # Root directory name
    root_dir = os.path.abspath(os.getcwd())
 
    # Directory with these properties
    dir_name = main_name + '_'.join(str(v) for v in properties.values())
    print('Processing ' + dir_name) 

    # Make the directory
    subprocess.run('mkdir ' + dir_name, shell=True, check=True)
    
    # Copy all the files from templates
    subprocess.run('cp -r templates/*.* ' + dir_name + '/', shell=True, check=True)
    
    # Go to directory
    os.chdir(dir_name)

    # Update parameters.txt
    sub_parameters(properties, param_file)
            
    # Compile 
    subprocess.run([py_version + ' ' + compile_script], shell=True)

    # Run in a screen       
    subprocess.call(['screen -d -m ./' + exe_name], shell=True)
                
    # Go back to root
    os.chdir(root_dir)

def run_simulations_serial(properties, param_file, main_name, py_version, compile_script, exe_name):
  
    # Root directory name
    root_dir = os.path.abspath(os.getcwd())
 
    # Directory with these properties
    dir_name = main_name + '_'.join(str(v) for v in properties.values())
    print('Processing ' + dir_name) 

    # Make the directory
    subprocess.run('mkdir ' + dir_name, shell=True, check=True)
    
    # Copy all the files from templates
    subprocess.run('cp -r templates/*.* ' + dir_name + '/', shell=True, check=True)
    
    # Go to directory
    os.chdir(dir_name)

    # Update parameters.txt
    sub_parameters(properties, param_file)
            
    # Compile 
    subprocess.run([py_version + ' ' + compile_script], shell=True)

    # Run        
    subprocess.call(['./' + exe_name], shell=True)
                
    # Go back to root
    os.chdir(root_dir)

def run_lj_simulation(
    sigma, epsilon, delta, r_cutoff, scale,
    param_file,      
    main_prefix,     
    py_version,      
    compile_script,  
    path_exe,        
    exe_name):
    """
    1) Update parameters.txt with LJ
    2) Compile C++ code
    3) Run the binary
    4) Rename results file to include r_cutoff
    """
    # Patch parameters
    sub_parameters({
        'Sigma':    sigma,
        'Epsilon':  epsilon,
        'Delta':    delta,
        'R_cutoff': r_cutoff,
        'Scale':    scale,
    }, param_file)

    # Compile
    subprocess.run([py_version, compile_script], check=True)

    # Run sim
    subprocess.run([os.path.join(path_exe, exe_name)], capture_output=True, check=True)

    # Rename output
    old = f"results_lj.txt"
    #suffix = f"{r_cutoff:.1e}".split("e")[0]
    suffix = f"{epsilon}_{delta}_{scale}"
    new = f"{main_prefix}{suffix}.txt"
    os.replace(old, new)

