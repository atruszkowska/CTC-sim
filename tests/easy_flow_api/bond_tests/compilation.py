# Script for compiling the tests

import subprocess, glob, os

### Input 
# Path to the main directory
path = '../../../src/'
# Path to executables 
path_exe = '../../../executables/'
# Compiler options
cx = 'g++'
std = '-std=c++11'
opt = '-O0'
other = '-Wall'
# Common source files
src_files = path + 'utils.cpp'
src_files += ' ' + path + 'test_utils.cpp'
src_files += ' ' + path + 'ctc_sim.cpp'
src_files += ' ' + path + 'boundaries/periodic_boundary.cpp'
src_files += ' ' + path + 'boundaries/wall_boundary.cpp'
src_files += ' ' + path + 'io_operations/file_handler.cpp'
src_files += ' ' + path + 'io_operations/load_parameters.cpp'
src_files += ' ' + path + 'interactions/basic.cpp'
src_files += ' ' + path + 'interactions/bond.cpp'
tst_files = ''

# Tests with bonds
# Name of the executable
exe_name = 'bnd_tst_api'
# Files needed only for this build
spec_files = 'api_bond_tests.cpp'
compile_command = ' '.join([cx, std, opt, other, '-o', exe_name, spec_files, tst_files, src_files])
subprocess.call([compile_command], shell=True)
subprocess.call(['mv ' + exe_name + ' ' + path_exe], shell=True)

# Tests without bonds
# Name of the executable
exe_name = 'nbnd_tst_api'
# Files needed only for this build
spec_files = 'api_bond_tests_no_bonds.cpp'
compile_command = ' '.join([cx, std, opt, other, '-o', exe_name, spec_files, tst_files, src_files])
subprocess.call([compile_command], shell=True)
subprocess.call(['mv ' + exe_name + ' ' + path_exe], shell=True)

# Bond tracking 
# Name of the executable
exe_name = 'bnd_tr_tst'
# Files needed only for this build
spec_files = 'bond_tracking_test.cpp'
compile_command = ' '.join([cx, std, opt, other, '-o', exe_name, spec_files, tst_files, src_files])
subprocess.call([compile_command], shell=True)
subprocess.call(['mv ' + exe_name + ' ' + path_exe], shell=True)


