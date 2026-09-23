# Test of cluster detection algorithms

import time, sys
import numpy as np
from collections import defaultdict
import math

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
t0 = '0'

# --- Data --------------------------------------------------------------------

pos_file = 'permanent_test_positions.txt'
bond_file = 'permanent_test_bonds.txt'

# --- Cluster detection -------------------------------------------------------    

df_pos = smut.read_pos(pos_file)
df_bond = smut.read_bonds(bond_file)

time0 = time.time()
uf_areas, uf_clusters, uf_counts = smut.cluster_data_union_find(df_bond, df_pos, Np, max_area, max_len, alpha)
time1 = time.time()

# Group data by realization number
realization_data = defaultdict(dict)
for key, value in uf_areas.items():
    time_step, ri = key.split('-')
    realization_data[ri][time_step] = value

realization_data_sizes = defaultdict(dict)
for key, value in uf_clusters.items():
    time_step, ri = key.split('-')
    realization_data_sizes[ri][time_step] = 0
    count = 0
    for outer in value:
        for cluster in outer:     
            realization_data_sizes[ri][time_step] += len(cluster)
            count += 1
    if count > 0:
        realization_data_sizes[ri][time_step] /= count
    else:
        realization_data_sizes[ri][time_step] = 0

realization_data_numbers = defaultdict(dict)
for key, value in uf_counts.items():
    time_step, ri = key.split('-')
    realization_data_numbers[ri][time_step] = value[0]

# Compute the averages
# Dictionary to store average ratios for each time step (other than t0)
average_ratios_by_time = {}
# Analogous dictionary for cluster sizes
average_cluster_sizes_by_time = {}
# Analogous dictionary for cluster numbers
average_cluster_numbers_by_time = {}

# Area ratios
# Loop through each realization
for ri, times in realization_data.items():  # ri is a string
    for tf, value in times.items():         # tf is a string
        if tf == t0:
            continue  # skip t0 itself
        if t0 in times and not math.isclose(float(times[t0]), 0):
            ratio = float(value) / float(times[t0])
            if ratio is not None and not math.isnan(ratio):
                if tf not in average_ratios_by_time:
                    average_ratios_by_time[tf] = []
                average_ratios_by_time[tf].append(ratio)

# Compute averages
for tf in average_ratios_by_time:
    ratios = average_ratios_by_time[tf]
    average_ratios_by_time[tf] = sum(ratios) / len(ratios) if ratios else float('nan') 
  

all_ts = {t for times in realization_data_sizes.values() for t in times}
average_cluster_sizes_by_time = {t: sum(times[t] for times in realization_data_sizes.values() if t in times) /
                      sum(1 for times in realization_data_sizes.values() if t in times) for t in all_ts}

all_ts = {t for times in realization_data_numbers.values() for t in times}
average_cluster_numbers_by_time = {t: sum(times[t] for times in realization_data_numbers.values() if t in times) /
                      sum(1 for times in realization_data_numbers.values() if t in times) for t in all_ts}

print(uf_clusters)
print(uf_counts)
print(uf_areas)
print(f"Time used {time1-time0}s")

print(realization_data)
print(realization_data_sizes)
print(realization_data_numbers)

print(average_ratios_by_time)
print(average_cluster_sizes_by_time)
print(average_cluster_numbers_by_time)
