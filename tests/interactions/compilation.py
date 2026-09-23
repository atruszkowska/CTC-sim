# Script for compiling the tests

import subprocess, glob, os

### Input 
# Path to the main directory
path = '../../src/'
# Path to executables 
path_exe = '../../executables/'
# Compiler options
cx = 'g++'
std = '-std=c++11'
opt = '-O0'
other = '-Wall'
# Common source files
src_files = path + 'utils.cpp'
src_files += ' ' + path + 'test_utils.cpp'
src_files += ' ' + path + 'io_operations/file_handler.cpp'
src_files += ' ' + path + 'interactions/basic.cpp'
src_files += ' ' + path + 'interactions/bond.cpp'
src_files += ' ' + path + 'interactions/wall_bond.cpp'
tst_files = ''

### Test suite 1
# Interactions tests
# Name of the executable
exe_name = 'int_test'
# Files needed only for this build
spec_files = 'interactions_tests.cpp'
compile_command = ' '.join([cx, std, opt, other, '-o', exe_name, spec_files, tst_files, src_files])
subprocess.call([compile_command], shell=True)
subprocess.call(['mv ' + exe_name + ' ' + path_exe], shell=True)


