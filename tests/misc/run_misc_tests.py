import subprocess

import sys
py_path = '../../scripts/'
sys.path.insert(0, py_path)

import utils as ut
from colors import *

py_version = 'python3'

# Directory with executables
path_exe = '../../executables/'

#
# Compile and run all the tests from misc category
#

# Compile
subprocess.call([py_version + ' compilation.py'], shell=True)

# Test utilities 
ut.msg('Test utility functions', CYAN)
subprocess.call([path_exe + 'tu_test'], shell=True)

# FileHandler class
ut.msg('File_handler class', CYAN)
subprocess.call([path_exe + 'file_hdl_tests'], shell=True)

# Easy_IO class
ut.msg('Easy_IO class', CYAN)
subprocess.call([path_exe + 'eio_tests'], shell=True)

# LoadParameters class
ut.msg('LoadParameters class', CYAN)
subprocess.call([path_exe + 'ldp_tests'], shell=True)
