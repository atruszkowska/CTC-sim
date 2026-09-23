import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

sys.path.insert(0, '../../scripts/')

import utils as ut
import simulation_utils as smut

# --- Matplotlib rcParams Update -----------------------------------------

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 14,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'xtick.major.size': 5,
    'ytick.major.size': 5,
    'svg.fonttype': 'none'
})

# --- System & Simulation Parameters ------------------------------------

Lx = 0.12
Ly = 1e-3 
max_area = Lx * Ly
max_len = 4.04e-5 
alpha = None
bond_gap = 2.2725e-05
t0 = '0'
delta_t = 0.1  # Converts step count to minutes

tau_values = [0, 2, 5, 20]
Np = 300

# Target plot snapshots in minutes
target_times_min = [0, 30, 60]

output_dir = 'key_cluster_snapshots'
os.makedirs(output_dir, exist_ok=True)

# Cluster color palette
cluster_colors = [
    '#e6194B', '#4363d8', '#911eb4', '#f032e6',
    '#800000', '#000075', '#469990', '#9A6324',
    '#e6beff', '#bfef45', '#f58231', '#3cb44b'
]

# --- Helper Functions ---------------------------------------------------

def get_df_step(df_pos, t_val, ri):
    """
    Safely extract position data for a specific timestep and realization index,
    matching column schemes ('it' and 'ir') or fallbacks ('timestep', 'realization').
    """
    cols_lower = {str(col).lower(): col for col in df_pos.columns}

    if 'it' in cols_lower and 'ir' in cols_lower:
        it_col, ir_col = cols_lower['it'], cols_lower['ir']
        df_step = df_pos[
            (df_pos[it_col] == float(t_val)) &
            (df_pos[ir_col].astype(str) == str(ri))
        ].copy()
    else:
        ts_col = next((cols_lower[c] for c in ['timestep', 'time_step', 'step', 't'] if c in cols_lower), None)
        ri_col = next((cols_lower[c] for c in ['realization', 'ri', 'realization_id', 'run'] if c in cols_lower), None)

        if ts_col and ri_col:
            df_step = df_pos[
                (df_pos[ts_col].astype(str) == str(t_val)) &
                (df_pos[ri_col].astype(str) == str(ri))
            ].copy()
        else:
            df_step = pd.DataFrame()

    if not df_step.empty:
        if 'it' in cols_lower and 'ir' in cols_lower:
            df_step['pid'] = df_step.groupby([cols_lower['it'], cols_lower['ir']]).cumcount()
        else:
            df_step['pid'] = range(len(df_step))
        df_step = df_step.sort_values('pid')

    return df_step


def generate_snapshot(df_step, clusters_for_step, title_str, save_path):
    """
    Plots monomers in grey (#d3d3d3) and multi-particle clusters in distinct colors
    with auto-scaling for optimal visual clarity.
    """
    if df_step is None or df_step.empty:
        return

    x_col, y_col = None, None
    for c in df_step.columns:
        c_str = str(c).strip().lower()
        if c_str in ['x', 'pos_x', 'x_pos', '0']:
            x_col = c
        elif c_str in ['y', 'pos_y', 'y_pos', '1']:
            y_col = c

    if x_col is None or y_col is None:
        numeric_cols = df_step.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) >= 2:
            x_col, y_col = numeric_cols[:2]
        else:
            return

    # Extract clusters of size > 1
    clusters = []
    for outer in clusters_for_step:
        for cluster in outer:
            if len(cluster) > 1:
                clusters.append(cluster)

    clustered_pids = {pid for cluster in clusters for pid in cluster}

    # Slim figure box for micro-channels
    fig, ax = plt.subplots(figsize=(12, 2.5))

    # Plot unclustered cells (monomers) in light grey
    monomers = [i for i in range(len(df_step)) if i not in clustered_pids]

    if monomers:
        ax.scatter(
            df_step.iloc[monomers][x_col],
            df_step.iloc[monomers][y_col],
            color='#d3d3d3',
            s=25,
            alpha=0.6,
            edgecolors='none',
            label='Monomers'
        )

    # Plot multi-particle clusters
    for i, cluster in enumerate(clusters):
        pids = [pid for pid in cluster if 0 <= pid < len(df_step)]
        if not pids:
            continue

        df_cl = df_step.iloc[pids]
        ax.scatter(
            df_cl[x_col],
            df_cl[y_col],
            color=cluster_colors[i % len(cluster_colors)],
            s=40,
            edgecolors='black',
            linewidth=0.5
        )

    # Bounds with auto aspect ratio to stretch y-axis comfortably
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    ax.set_aspect('auto')

    ax.set_xlabel('X Position (m)')
    ax.set_ylabel('Y Position (m)')
    ax.set_title(title_str, pad=10)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close(fig)


# --- Main Processing Loop -----------------------------------------------

for tau in tau_values:

    pos_file = f'results{tau}.txt'
    bond_file = f'bonds{tau}.txt'

    if not os.path.exists(pos_file) or not os.path.exists(bond_file):
        print(f"Skipping TAU = {tau}: file(s) {pos_file} or {bond_file} not found.")
        continue

    # 1. Ingest dataset
    df_pos = smut.read_pos(pos_file)
    df_bond = smut.read_bonds(bond_file)
    
    # 2. Extract cluster topology using Union-Find
    uf_areas, uf_clusters, uf_counts = smut.cluster_data_union_find(
        df_bond, df_pos, Np, max_area, max_len, alpha
    )

    # 3. Group area values by realization number
    realization_data = defaultdict(dict)
    for key, value in uf_areas.items():
        time_step, ri = key.split('-')
        realization_data[ri][time_step] = value

    # 4. Calculate final area ratios (A_final / A0)
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
        print(f"No valid realizations found for TAU = {tau}")
        continue

    # 5. Identify the 4 Key Realizations
    sorted_ris = sorted(final_ratios.keys(), key=lambda r: final_ratios[r])
    mean_ratio = np.mean(list(final_ratios.values()))

    min_ri = sorted_ris[0]
    max_ri = sorted_ris[-1]
    target_ri = sorted_ris[len(sorted_ris) // 2]
    avg_ri = min(final_ratios.keys(), key=lambda r: abs(final_ratios[r] - mean_ratio))

    key_realizations = {
        'min': min_ri,
        'max': max_ri,
        'target': target_ri,
        'avg': avg_ri
    }

    print(f"\n--- Key Realizations for TAU = {tau} ---")
    for cat, r_id in key_realizations.items():
        print(f"  {cat.capitalize():<6}: Realization {r_id} (Final Area Ratio: {final_ratios[r_id]:.4f})")

    # 6. Generate snapshots at target times (0, 30, 60 min)
    for cat_name, ri in key_realizations.items():
        if ri not in realization_data:
            continue

        available_steps = sorted(realization_data[ri].keys(), key=float)
        available_mins = np.array([float(s) * delta_t / 60.0 for s in available_steps])

        for target_m in target_times_min:
            closest_idx = np.argmin(np.abs(available_mins - target_m))
            t_val = available_steps[closest_idx]
            actual_m = available_mins[closest_idx]

            t_float = float(t_val)
            formatted_t = f"{int(t_float)}" if t_float.is_integer() else f"{t_val}"
            
            lookup_key = f"{t_val}-{ri}"
            clusters_for_step = uf_clusters.get(lookup_key, uf_clusters.get(f"{formatted_t}-{ri}", []))

            df_step = get_df_step(df_pos, t_val, ri)
            if df_step.empty:
                continue

            title = f'TAU {tau} | Case: {cat_name.upper()} (Realization {ri}) | Time: {target_m} min ({actual_m:.1f} min)'
            out_filename = os.path.join(
                output_dir,
                f'snapshot_tau_{tau}_{cat_name}_ri_{ri}_t_{target_m}min.png'
            )

            generate_snapshot(df_step, clusters_for_step, title, out_filename)

print("\nSnapshots rendered with clear auto-scaling!")
