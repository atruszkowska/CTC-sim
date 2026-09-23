import os
import sys
import math
import re
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib.ticker import MaxNLocator

sys.path.insert(0, '../../scripts/')

import utils as ut
import simulation_utils as smut
import visualization as vis

# --- Matplotlib rcParams Update -----------------------------------------

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 18,
    'axes.labelsize': 18,
    'axes.titlesize': 18,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 14,
    'axes.linewidth': 1.5,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'xtick.major.size': 6,
    'ytick.major.size': 6,
    'svg.fonttype': 'none'
})

# --- Input parameters ---------------------------------------------------

Lx = 0.12
Ly = 1e-3 
max_area = Lx * Ly
max_len = 4.04e-5 
alpha = None
bond_gap = 2.2725e-05
t0 = '0'
delta_t = 0.1

tau_values = [0, 2, 5, 20]
data_dir = '.'

y_min = 0.0
y_max = {'Cluster size': 4, 'Number of clusters': 75, r'Area ratio ($A(t)/A_0$)': 2}

# Variables configurations (color, y-label, file suffix)
variables_config = {
    'ratio': {
        'color': 'darkred',
        'ylabel': r'Area ratio ($A(t)/A_0$)',
        'filename_suffix': 'area_ratio'
    },
    'size': {
        'color': 'purple',
        'ylabel': 'Cluster size',
        'filename_suffix': 'cluster_size'
    },
    'number': {
        'color': 'navy',
        'ylabel': 'Number of clusters',
        'filename_suffix': 'cluster_number'
    }
}

# Line styles for key realizations
linestyles = {
    'target': ('-', 2.5),
    'min': ('--', 2.0),
    'max': (':', 2.5),
    'avg': ('-.', 2.0)
}

# --- Find Valid File Pairs (Base & Number_of_particles_200/400) ---------

results_files = {}

# Matches strictly:
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

# --- Main Processing Loop -----------------------------------------------

for (wss, suffix), (pos_file, Np) in sorted(results_files.items()):

    bond_file = f'bonds{wss}_{suffix}.txt'
    pos_path = os.path.join(data_dir, pos_file)
    bond_path = os.path.join(data_dir, bond_file)

    if not os.path.exists(bond_path):
        print(f"Skipping WSS = {wss} | {suffix}: missing file {bond_file}")
        continue

    print(f"Processing WSS = {wss} | {suffix} (Np = {Np})")

    df_pos = smut.read_pos(pos_path)
    df_bond = smut.read_bonds(bond_path)
    
    uf_areas, uf_clusters, uf_counts = smut.cluster_data_union_find(
        df_bond, df_pos, Np, max_area, max_len, alpha
    )
 
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

    # Calculate final area ratio for each valid realization
    final_ratios = {}
    for ri, times in realization_data.items():
        if t0 not in times:
            continue

        t_last = max(times.keys(), key=float)
        A0 = float(times[t0])

        if math.isclose(A0, 0):
            continue

        final_ratios[ri] = float(times[t_last]) / A0

    if not final_ratios:
        print(f"No valid realizations found for WSS = {wss} | {suffix}")
        continue

    # --- Build time series for all realizations ----------------------------

    time_series = {}

    for ri in realization_data:
        if t0 not in realization_data[ri]:
            continue

        t_steps = sorted(realization_data[ri].keys(), key=float)
        times_float = [float(t) * delta_t / 60 for t in t_steps]

        A0 = float(realization_data[ri][t0])
        if math.isclose(A0, 0):
            continue

        time_series[ri] = {
            'times': times_float,
            'ratio': [float(realization_data[ri][t]) / A0 for t in t_steps],
            'size': [realization_data_sizes[ri][t] for t in t_steps],
            'number': [realization_data_numbers[ri][t] for t in t_steps]
        }

    if not time_series:
        print(f"No valid realizations found for WSS = {wss} | {suffix}")
        continue

    # Assume all realizations have the same time steps
    first_ri = next(iter(time_series))
    times = np.array(time_series[first_ri]['times'])

    # Calculate ensemble averages
    avg_series = {}
    for var_key in ['ratio', 'size', 'number']:
        values = np.array([
            time_series[ri][var_key]
            for ri in time_series
        ])
        avg_series[var_key] = np.mean(values, axis=0)

    # Save ensemble averages
    with open(f'avg_numbers_{wss}_{suffix}.txt', 'w') as f:
        for t, value in zip(times, avg_series['number']):
            f.write(f'{t} {value}\n')

    with open(f'avg_sizes_{wss}_{suffix}.txt', 'w') as f:
        for t, value in zip(times, avg_series['size']):
            f.write(f'{t} {value}\n')

    with open(f'avg_ratios_{wss}_{suffix}.txt', 'w') as f:
        for t, value in zip(times, avg_series['ratio']):
            f.write(f'{t} {value}\n')

    # --- Plotting: all realizations + ensemble average --------------------

    for var_key, config in variables_config.items():

        fig, ax = plt.subplots(figsize=(7.5, 5.5))

        color = config['color']
        ylabel = config['ylabel']

        # Plot individual realizations
        for ri, data in time_series.items():
            ax.plot(
                data['times'], data[var_key],
                color=color,
                alpha=0.08,
                linewidth=0.8
            )

        # Plot ensemble average
        ax.plot(
            times,
            avg_series[var_key],
            color=color,
            linewidth=3.0,
            label='Average'
        )

        ax.set_xlabel('Time (min)')
        ax.set_ylabel(ylabel)

        ax.set_xlim(0, 60)
        ax.set_ylim(y_min, y_max[ylabel])

        locator = MaxNLocator(nbins=6)
        ax.yaxis.set_major_locator(locator)

        ax.legend(loc='best', framealpha=0.9)
        plt.tight_layout()

        filename = f'tau_{wss}_{suffix}_{config["filename_suffix"]}_time_evolution.svg'
        plt.savefig(filename)
        plt.close()
