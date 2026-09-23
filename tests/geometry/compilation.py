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
tst_files = ''

### Test suite 1
# Box class 
# Name of the executable
exe_name = 'bx_test'
# Files needed only for this build
spec_files = 'box_tests.cpp'
compile_command = ' '.join([cx, std, opt, other, '-o', exe_name, spec_files, tst_files, src_files])
subprocess.call([compile_command], shell=True)
subprocess.call(['mv ' + exe_name + ' ' + path_exe], shell=True)

### Test suite 2
# Geometry class 
# Name of the executable
exe_name = 'gm_test'
# Files needed only for this build
spec_files = 'geometry_tests.cpp'
compile_command = ' '.join([cx, std, opt, other, '-o', exe_name, spec_files, tst_files, src_files])
subprocess.call([compile_command], shell=True)
subprocess.call(['mv ' + exe_name + ' ' + path_exe], shell=True)

### Test suite 3
# Circle class 
# Name of the executable
exe_name = 'cr_test'
# Files needed only for this build
spec_files = 'circle_tests.cpp'
compile_command = ' '.join([cx, std, opt, other, '-o', exe_name, spec_files, tst_files, src_files])
subprocess.call([compile_command], shell=True)
subprocess.call(['mv ' + exe_name + ' ' + path_exe], shell=True)


