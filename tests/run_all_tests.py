import subprocess, os

import sys
py_path = '../scripts/'
sys.path.insert(0, py_path)

import utils as ut
from colors import *

py_version = 'python3'

# 
# Compile and run all the tests
#

# Bar delimiting the test suites
nSim = 10

print('\n'*2)

# Misc
print('\n'*2)
ut.msg('- '*nSim + 'MISCELLANEOUS TESTS' + ' -'*nSim, REVERSE+RED)
os.chdir('misc/')
subprocess.call([py_version + ' run_misc_tests.py'], shell=True)
os.chdir('../')

# Random number generator class
print('\n'*2)
ut.msg('- '*nSim + 'RNG CLASS TESTS' + ' -'*nSim, REVERSE+RED)
os.chdir('rng/')
subprocess.call([py_version + ' run_rng_tests.py'], shell=True)
os.chdir('../')

# Easy Flow Api
print('\n'*2)
ut.msg('- '*nSim + 'EASY FLOW API TESTS' + ' -'*nSim, REVERSE+RED)
os.chdir('ctc_sim_api/')
subprocess.call([py_version + ' run_api_tests.py'], shell=True)
os.chdir('../')

# Boundaries 
print('\n'*2)
ut.msg('- '*nSim + 'BOUNDARIES TESTS' + ' -'*nSim, REVERSE+RED)
os.chdir('boundaries/')
subprocess.call([py_version + ' run_boundaries_tests.py'], shell=True)
os.chdir('../')

# Geometry class
print('\n'*2)
ut.msg('- '*nSim + 'GEOMETRY TESTS' + ' -'*nSim, REVERSE+RED)
os.chdir('geometry/')
subprocess.call([py_version + ' run_geometry_tests.py'], shell=True)
os.chdir('../')

# Interactions clss
print('\n'*2)
ut.msg('- '*nSim + 'INTERACTIONS TESTS' + ' -'*nSim, REVERSE+RED)
os.chdir('interactions/')
subprocess.call([py_version + ' run_interactions_tests.py'], shell=True)
os.chdir('../')

# Particle class
print('\n'*2)
ut.msg('- '*nSim + 'PARTICLES TESTS' + ' -'*nSim, REVERSE+RED)
os.chdir('particles/')
subprocess.call([py_version + ' run_particle_tests.py'], shell=True)
os.chdir('../')
