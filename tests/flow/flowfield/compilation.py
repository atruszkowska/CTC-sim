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
src_files = path + 'flow.cpp' + ' ' + path + 'utils.cpp '
src_files += ' ' + path + 'test_utils.cpp'
tst_files = ''

# Name of the executable
exe_name = 'flf_test'
# Files needed only for this build
spec_files = 'flow_field_tests.cpp '
compile_com = ' '.join([cx, std, opt, other, '-o', exe_name, spec_files, tst_files, src_files])
subprocess.call([compile_com], shell=True)
subprocess.call(['mv ' + exe_name + ' ' + path_exe], shell=True)


