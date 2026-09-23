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
ut.msg('Boundary conditions tests', CYAN)
subprocess.call([path_exe + 'bc_test'], shell=True)

ut.msg('Boundary xy tests', CYAN)
subprocess.call([path_exe + 'bcxy_test'], shell=True)
