import subprocess

import sys
py_path = '../../../scripts/'
sys.path.insert(0, py_path)

import utils as ut
from colors import *

py_version = 'python3'

# Directory with executables
path_exe = '../../../executables/'

#
# Compile and run all the FlowField class specific tests
#

# Compile
subprocess.call([py_version + ' compilation.py'], shell=True)

# Run
ut.msg('FlowField class tests', CYAN)
subprocess.call([path_exe + 'flfunc_test'], shell=True)


