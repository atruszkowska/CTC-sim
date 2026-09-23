# Script for compiling the simulations

import subprocess, glob, os

### Input 
# Path to the main directory
path = '../../src/'
# Path to executables 
path_exe = '../../executables/'
# Compiler options
cx = 'g++'
std = '-std=c++11'
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
src_files += ' ' + '../../include/' + 'io_operations/easy_io.h'
src_files += ' ' + path + 'io_operations/load_parameters.cpp'
src_files += ' ' + path + 'interactions/basic.cpp'
src_files += ' ' + path + 'interactions/bond.cpp'
tst_files = ''

# Name of the executable
exe_name = 'bc_test'
# Files needed only for this build
spec_files = 'bc_tests.cpp'
compile_command = ' '.join([cx, std, opt, other, '-o', exe_name, spec_files, tst_files, src_files])
subprocess.call([compile_command], shell=True)
subprocess.call(['mv ' + exe_name + ' ' + path_exe], shell=True)

# Name of the executable
exe_name = 'bcxy_test'
# Files needed only for this build
spec_files = 'bc_xy_tests.cpp'
compile_command = ' '.join([cx, std, opt, other, '-o', exe_name, spec_files, tst_files, src_files])
subprocess.call([compile_command], shell=True)
subprocess.call(['mv ' + exe_name + ' ' + path_exe], shell=True)
