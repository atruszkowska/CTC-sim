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
src_files += ' ' + path + 'io_operations/load_parameters.cpp'
tst_files = ''

### Test suite 1
# test_utils tests
# Name of the executable
exe_name = 'tu_test'
# Files needed only for this build
spec_files = 'test_utils_tests.cpp'
compile_command = ' '.join([cx, std, opt, other, '-o', exe_name, spec_files, tst_files, src_files])
subprocess.call([compile_command], shell=True)
subprocess.call(['mv ' + exe_name + ' ' + path_exe], shell=True)

### Test suite 2
# File_handler.h tests
# ------ Extra condition 
# Remove the custom mode test file if exists
f_app = './test_data/custom_mode.txt'
if os.path.exists(f_app): 
	os.remove(f_app)
# Name of the executable
exe_name = 'file_hdl_tests'
# Files needed only for this build
spec_files = 'file_handler_tests.cpp'
compile_com = ' '.join([cx, std, opt, '-o', exe_name, spec_files, src_files, tst_files])
subprocess.call([compile_com], shell=True)
subprocess.call(['mv ' + exe_name + ' ' + path_exe], shell=True)

# easy_io.h tests
# Remove these files if exist
f_path = './test_data/'
files_rm = ['wr_bool.txt', 'wr_int.txt', 'wr_string.txt', 'wr_double.txt', 'wr_bool_1D.txt', 'wr_int_1D.txt', 'wr_string_1D.txt', 'wr_double_1D.txt']
files_rm = [f_path + x for x in files_rm]
files_rm.append('./dflt_abm_io_file.txt')
for frm in files_rm:
	if os.path.exists(frm): 
		os.remove(frm)
# Name of the executable
exe_name = 'eio_tests'
# Files needed only for this build
spec_files = 'easy_io_tests.cpp'
compile_com = ' '.join([cx, std, opt, '-o', exe_name, src_files, spec_files, tst_files])
subprocess.call([compile_com], shell=True)
subprocess.call(['mv ' + exe_name + ' ' + path_exe], shell=True)

# load_parameter.h tests
# Remove these files if exist
f_path = './test_data/'
files_rm = ['wr_bool.txt', 'wr_int.txt', 'wr_string.txt', 'wr_double.txt', 'wr_bool_1D.txt', 'wr_int_1D.txt', 'wr_string_1D.txt', 'wr_double_1D.txt']
files_rm = [f_path + x for x in files_rm]
files_rm.append('./dflt_abm_io_file.txt')
for frm in files_rm:
	if os.path.exists(frm): 
		os.remove(frm)
# Name of the executable
exe_name = 'ldp_tests'
# Files needed only for this build
spec_files = 'load_parameters_tests.cpp'
compile_com = ' '.join([cx, std, opt, '-o', exe_name, spec_files, src_files, tst_files])
subprocess.call([compile_com], shell=True)
subprocess.call(['mv ' + exe_name + ' ' + path_exe], shell=True)


