import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib.ticker import MaxNLocator

sys.path.insert(0, '../../scripts/')

import utils as ut
import simulation_utils as smut
import visualization as vis

# ============================================================
# USER SETTINGS
# ============================================================

tau_values = [0, 2, 5, 20]

t0 = '0'
delta_t = 0.1
t1 = int(15 * 60 / delta_t)       # 15 min = 9000 time steps

Np = 300
Lx = 0.12
Ly = 1e-3
max_area = Lx * Ly
max_len = 4.04e-5
alpha = None

output_file = 'temporal_statistics_all_wss.txt'

# ============================================================
# OUTPUT FILE
# ============================================================

with open(output_file, 'w') as fout:

    fout.write(
        'Temporal averages from 15 min to end\n'
        'Area ratio includes all valid time steps.\n'
        'Number of clusters and cluster size include only time steps '
        'where clusters exist (N > 0).\n\n'
    )

    # ========================================================
    # LOOP OVER WALL SHEAR STRESSES
    # ========================================================

    for wss in tau_values:

        print(f'\nWSS = {wss} dyn/cm^2')

        # ----------------------------------------------------
        # READ DATA
        # ----------------------------------------------------

        df_pos = smut.read_pos(f'results{wss}.txt')
        df_bond = smut.read_bonds(f'bonds{wss}.txt')

        uf_areas, uf_clusters, uf_counts = smut.cluster_data_union_find(
            df_bond,
            df_pos,
            Np,
            max_area,
            max_len,
            alpha
        )

        # ----------------------------------------------------
        # GROUP DATA BY REALIZATION
        # ----------------------------------------------------

        realization_areas = defaultdict(dict)

        for key, value in uf_areas.items():
            time_step, ri = key.split('-')
            realization_areas[ri][time_step] = value

        realization_sizes = defaultdict(dict)

        for key, value in uf_clusters.items():

            time_step, ri = key.split('-')

            realization_sizes[ri][time_step] = 0
            count = 0

            for outer in value:
                for cluster in outer:
                    realization_sizes[ri][time_step] += len(cluster)
                    count += 1

            if count > 0:
                realization_sizes[ri][time_step] /= count
            else:
                realization_sizes[ri][time_step] = 0

        realization_numbers = defaultdict(dict)

        for key, value in uf_counts.items():
            time_step, ri = key.split('-')
            realization_numbers[ri][time_step] = value[0]

        # ----------------------------------------------------
        # CALCULATE TEMPORAL AVERAGE FOR EACH REALIZATION
        # ----------------------------------------------------

        realization_stats = {}

        for ri in realization_areas:

            # Need the initial area
            if t0 not in realization_areas[ri]:
                continue

            A0 = float(realization_areas[ri][t0])

            if math.isclose(A0, 0):
                continue

            area_values = []
            size_values = []
            number_values = []

            # -----------------------------------------------
            # Loop over time steps
            # -----------------------------------------------

            for t in realization_areas[ri]:

                if int(t) < t1:
                    continue

                area = realization_areas[ri][t]

                if isinstance(area, float) and math.isnan(area):
                    continue

                # Area ratio:
                # ALL valid time steps after 15 min
                area_values.append(float(area) / A0)

                # Number and size:
                # ONLY time steps where clusters exist
                if t in realization_numbers[ri]:

                    number = realization_numbers[ri][t]

                    if number > 0:

                        number_values.append(float(number))

                        if t in realization_sizes[ri]:
                            size_values.append(
                                float(realization_sizes[ri][t])
                            )

            # No valid area data
            if not area_values:
                continue

            realization_stats[ri] = {
                'ratio': np.mean(area_values),
                'size': (
                    np.mean(size_values)
                    if size_values
                    else np.nan
                ),
                'number': (
                    np.mean(number_values)
                    if number_values
                    else np.nan
                )
            }

        # ----------------------------------------------------
        # ENSEMBLE STATISTICS
        # ----------------------------------------------------

        ratio_values = []
        number_values = []
        size_values = []

        for ri, stats in realization_stats.items():

            ratio_values.append(stats['ratio'])

            if not np.isnan(stats['number']):
                number_values.append(stats['number'])

            if not np.isnan(stats['size']):
                size_values.append(stats['size'])

        # ----------------------------------------------------
        # CALCULATE MEAN / STD / MIN / MAX
        # ----------------------------------------------------

        if ratio_values:
            ratio_mean = np.mean(ratio_values)
            ratio_std = np.std(ratio_values)
            ratio_min = np.min(ratio_values)
            ratio_max = np.max(ratio_values)
        else:
            ratio_mean = ratio_std = ratio_min = ratio_max = np.nan

        if number_values:
            number_mean = np.mean(number_values)
            number_std = np.std(number_values)
            number_min = np.min(number_values)
            number_max = np.max(number_values)
        else:
            number_mean = number_std = number_min = number_max = np.nan

        if size_values:
            size_mean = np.mean(size_values)
            size_std = np.std(size_values)
            size_min = np.min(size_values)
            size_max = np.max(size_values)
        else:
            size_mean = size_std = size_min = size_max = np.nan

        # ----------------------------------------------------
        # PRINT
        # ----------------------------------------------------

        print(
            f'  Area ratio:       '
            f'{ratio_mean:.6g} +/- {ratio_std:.6g} '
            f'[{ratio_min:.6g}, {ratio_max:.6g}]'
        )

        print(
            f'  Number clusters:  '
            f'{number_mean:.6g} +/- {number_std:.6g} '
            f'[{number_min:.6g}, {number_max:.6g}]'
        )

        print(
            f'  Cluster size:     '
            f'{size_mean:.6g} +/- {size_std:.6g} '
            f'[{size_min:.6g}, {size_max:.6g}]'
        )

        print(
            f'  Realizations:     {len(realization_stats)}'
        )

        print(
            f'  N realizations:   {len(number_values)}'
        )

        print(
            f'  S realizations:   {len(size_values)}'
        )

        # ----------------------------------------------------
        # WRITE TO FILE
        # ----------------------------------------------------

        fout.write(f'WSS = {wss} dyn/cm^2\n')

        fout.write(
            f'Area ratio: '
            f'mean = {ratio_mean:.6g}, '
            f'std = {ratio_std:.6g}, '
            f'min = {ratio_min:.6g}, '
            f'max = {ratio_max:.6g}\n'
        )

        fout.write(
            f'Number of clusters: '
            f'mean = {number_mean:.6g}, '
            f'std = {number_std:.6g}, '
            f'min = {number_min:.6g}, '
            f'max = {number_max:.6g}\n'
        )

        fout.write(
            f'Cluster size: '
            f'mean = {size_mean:.6g}, '
            f'std = {size_std:.6g}, '
            f'min = {size_min:.6g}, '
            f'max = {size_max:.6g}\n'
        )

        fout.write(
            f'Number of realizations with area data: '
            f'{len(realization_stats)}\n'
        )

        fout.write(
            f'Number of realizations with clusters: '
            f'{len(number_values)}\n'
        )

        fout.write(
            f'Number of realizations with cluster-size data: '
            f'{len(size_values)}\n'
        )

        fout.write('\n')

print(f'\nResults written to {output_file}')
