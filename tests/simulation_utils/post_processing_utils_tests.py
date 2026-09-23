from pathlib import Path
import sys

import pandas as pd
sys.path.insert(0, '../../scripts/')
import utils as ut
import simulation_utils as smut
from colors import *

num_particles = 4
time_steps = [0, 1]
reals = [0, 1]

pos_rows = []
for ir in reals:
    for it in time_steps:
        for pid in range(num_particles):
            y_coord = pid & 1
            pos_rows.append((pid, y_coord, it, ir))
pos_df = pd.DataFrame(pos_rows, columns=["x", "y", "it", "ir"])

bond_rows = [
    (0, 1, 0, 0), (2, 3, 0, 0),
    (0, 1, 1, 0),
    (0, 1, 0, 1), (1, 2, 0, 1),
    (0, 1, 1, 1), (1, 2, 1, 1), (2, 3, 1, 1)
]
bond_df = pd.DataFrame(bond_rows, columns=["ID-1", "ID-2", "it", "ir"])

max_area = 10.0
max_len = 10.0
mean_areas, cluster_lists, cluster_counts = smut.cluster_data_union_find(
    bond_df, pos_df, num_particles, max_area, max_len, alpha=0.1
)

expected_clusters = {
    "0-0": [[[0, 1], [2, 3]]],
    "1-0": [[[0, 1]]],
    "0-1": [[[0, 1, 2]]],
    "1-1": [[[0, 1, 2, 3]]],
}
expected_counts = {
    "0-0": [2],
    "1-0": [1],
    "0-1": [1],
    "1-1": [1],
}
sqrt2 = 2 ** 0.5
expected_areas = {
    "0-0": sqrt2,
    "1-0": sqrt2,
    "0-1": 1.0,
    "1-1": 2.0,
}

ok_lists = pd.Series(cluster_lists, dtype=object).sort_index().equals(
    pd.Series(expected_clusters, dtype=object).sort_index()
)
ok_counts = pd.Series(cluster_counts, dtype=object).sort_index().equals(
    pd.Series(expected_counts, dtype=object).sort_index()
)
ok_areas = pd.Series(mean_areas).sort_index().round(3).equals(
    pd.Series(expected_areas).sort_index().round(3)
)

ut.test_pass(ok_lists, "Cluster lists")
ut.test_pass(ok_counts, "Cluster counts")
ut.test_pass(ok_areas, "Cluster areas")

