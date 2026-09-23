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
# Compile and run the bond test suite 
#

# Compile
subprocess.call([py_version + ' compilation.py'], shell=True)

# Test utilities 
ut.msg('Cell bonds - API', CYAN)
subprocess.call([path_exe + 'bnd_tst_api'], shell=True)

ut.msg('Cell bonds - API - no bonds', CYAN)
subprocess.call([path_exe + 'nbnd_tst_api'], shell=True)

ut.msg('API tracking of cell bonds', CYAN)
subprocess.call([path_exe + 'bnd_tr_tst'], shell=True)


