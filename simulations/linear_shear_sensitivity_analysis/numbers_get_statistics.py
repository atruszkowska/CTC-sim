import os
import sys
import math
import re
import numpy as np
from collections import defaultdict

sys.path.insert(0, '../../scripts/')

import simulation_utils as smut


# ============================================================
# USER SETTINGS
# ============================================================

tau_values = [0, 2, 5, 20]

t0 = '0'
delta_t = 0.1
t1 = int(15 * 60 / delta_t)       # 15 min = 9000 time steps

Lx = 0.12
Ly = 1e-3
max_area = Lx * Ly
max_len = 4.04e-5
alpha = None

# Directory containing the results/bonds files
data_dir = '.'

output_file = 'temporal_statistics_particle_counts.txt'


# ============================================================
# FIND FILE PAIRS (Base & Number_of_particles_200/400)
# ============================================================

results_files = {}

# Matches:
# 1. results<wss>_base.txt
# 2. results<wss>_Number_of_particles_<200|400>.txt
pattern = re.compile(
    r'^results(0|2|5|20)_(?:(base)|Number_of_particles_(200|400))\.txt$'
)

for filename in os.listdir(data_dir):
    match = pattern.match(filename)
    if not match:
        continue

    wss = int(match.group(1))
    is_base = match.group(2)
    particles = match.group(3)

    if is_base:
        suffix = 'base'
        np_val = 300
    else:
        suffix = f'Number_of_particles_{particles}'
        np_val = int(particles)

    results_files[(wss, suffix)] = (filename, np_val)


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
    # LOOP OVER MATCHED FILES
    # ========================================================

    for (wss, suffix), (results_file, Np) in sorted(results_files.items()):

        bonds_file = f'bonds{wss}_{suffix}.txt'
        bonds_path = os.path.join(data_dir, bonds_file)

        if not os.path.exists(bonds_path):
            print(
                f'WARNING: missing bond file for '
                f'{results_file}: {bonds_file}'
            )
            continue

        print(
            f'\nWSS = {wss} dyn/cm^2'
            f' | {suffix} (Np = {Np})'
        )

        # ----------------------------------------------------
        # READ DATA
        # ----------------------------------------------------

        results_path = os.path.join(data_dir, results_file)

        df_pos = smut.read_pos(results_path)
        df_bond = smut.read_bonds(bonds_path)

        uf_areas, uf_clusters, uf_counts = (
            smut.cluster_data_union_find(
                df_bond,
                df_pos,
                Np,
                max_area,
                max_len,
                alpha
            )
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

            if t0 not in realization_areas[ri]:
                continue

            A0 = float(realization_areas[ri][t0])

            if math.isclose(A0, 0):
                continue

            area_values = []
            size_values = []
            number_values = []

            for t in realization_areas[ri]:

                if int(t) < t1:
                    continue

                area = realization_areas[ri][t]

                if isinstance(area, float) and math.isnan(area):
                    continue

                area_values.append(float(area) / A0)

                if t in realization_numbers[ri]:
                    number = realization_numbers[ri][t]

                    if number > 0:
                        number_values.append(float(number))

                        if t in realization_sizes[ri]:
                            size_values.append(
                                float(realization_sizes[ri][t])
                            )

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
        # PRINT & WRITE TO FILE
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

        fout.write(f'WSS = {wss} dyn/cm^2\n')
        fout.write(f'Parameter set = {suffix} (Np = {Np})\n')
        fout.write(
            f'Area ratio: mean = {ratio_mean:.6g}, std = {ratio_std:.6g}, '
            f'min = {ratio_min:.6g}, max = {ratio_max:.6g}\n'
        )
        fout.write(
            f'Number of clusters: mean = {number_mean:.6g}, std = {number_std:.6g}, '
            f'min = {number_min:.6g}, max = {number_max:.6g}\n'
        )
        fout.write(
            f'Cluster size: mean = {size_mean:.6g}, std = {size_std:.6g}, '
            f'min = {size_min:.6g}, max = {size_max:.6g}\n'
        )
        fout.write(
            f'Number of realizations with area data: {len(realization_stats)}\n'
        )
        fout.write(
            f'Number of realizations with clusters: {len(number_values)}\n'
        )
        fout.write(
            f'Number of realizations with cluster-size data: {len(size_values)}\n\n'
        )


print(f'\nResults written to {output_file}')
