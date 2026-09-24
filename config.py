#
# Sets up most of the components needed to use this software
#

import subprocess, os
from pathlib import Path

# Make a directory to store executables (if it doesn't already exists)
exe_dir = 'executables'
if not os.path.isdir(exe_dir):
	subprocess.call(['mkdir ' + exe_dir], shell=True)

# Make a directory to store to store documentation 
doc_dir = 'documentation'
if not os.path.isdir(doc_dir):
	subprocess.call(['mkdir ' + doc_dir], shell=True)


