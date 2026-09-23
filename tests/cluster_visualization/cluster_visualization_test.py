# Test of cluster visualization functions

import time, sys
import numpy as np

sys.path.insert(0, '../../scripts/')
import utils as ut
import simulation_utils as smut
import visualization as vis
 
# --- Input parameters ---------------------------------------------------

Np = 100 
Lx = 0.12
Ly = 1e-3 
max_area = Lx*Ly
max_len = np.hypot(Lx, Ly)
alpha = None
bond_gap = 2.02e-5
particle_radius = 1.01e-5
times_to_plot = [0, 6, 10]
rep_ir = '0'
# --- Data --------------------------------------------------------------------

pos_file = '../cluster_detection/permanent_test_positions.txt'
bond_file = '../cluster_detection/permanent_test_bonds.txt'

# --- Cluster detection -------------------------------------------------------    
import pandas as pd

df_pos = smut.read_pos(pos_file)
df_bond = smut.read_bonds(bond_file)

uf_areas, uf_clusters, uf_counts = smut.cluster_data_union_find(df_bond, df_pos, Np, max_area, max_len, alpha)

vis.plot_cells_clusters(df_pos, uf_clusters, rep_ir, times_to_plot, Lx, Ly, particle_radius)
