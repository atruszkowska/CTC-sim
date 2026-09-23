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
# Compile and run all the Flow class specific tests
#

# Compile
subprocess.call([py_version + ' compilation.py'], shell=True)

# Run
ut.msg('Force class tests', CYAN)
subprocess.call([path_exe + 'fc_test'], shell=True)


