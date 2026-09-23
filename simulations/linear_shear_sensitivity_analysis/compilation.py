# Script for compiling the simulations

import subprocess, glob, os

# exe_name is a commandline argument

### Input 
# Path to the main directory
path = '../../src/'
# Compiler options
cx = 'g++'
std = '-std=c++17'
opt = '-O3'
other = '-Wall'
# Common source files
src_files = path + 'utils.cpp'
src_files += ' ' + path + 'test_utils.cpp'
src_files += ' ' + path + 'ctc_sim.cpp'
src_files += ' ' + path + 'boundaries/periodic_boundary.cpp'
src_files += ' ' + path + 'boundaries/wall_boundary.cpp'
src_files += ' ' + path + 'boundaries/lj_wall.cpp'
src_files += ' ' + path + 'boundaries/round_wall.cpp'
src_files += ' ' + path + 'io_operations/file_handler.cpp'
src_files += ' ' + path + 'io_operations/load_parameters.cpp'
src_files += ' ' + path + 'interactions/basic.cpp'
src_files += ' ' + path + 'interactions/bond.cpp'
src_files += ' ' + path + 'interactions/linear_elastic_bond.cpp'
src_files += ' ' + path + 'interactions/wall_bond.cpp'
tst_files = ''

# Files needed only for this build
spec_files = 'simulation.cpp'
exe_name = 'lns'
compile_command = ' '.join([cx, std, opt, other, '-o', exe_name, spec_files, tst_files, src_files])
subprocess.call([compile_command], shell=True)


