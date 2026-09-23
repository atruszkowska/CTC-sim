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

# Run tests  
ut.msg('Box tests', CYAN)
subprocess.call([path_exe + 'bx_test'], shell=True)

ut.msg('Circle tests', CYAN)
subprocess.call([path_exe + 'cr_test'], shell=True)

ut.msg('Geometry tests', CYAN)
subprocess.call([path_exe + 'gm_test'], shell=True)

