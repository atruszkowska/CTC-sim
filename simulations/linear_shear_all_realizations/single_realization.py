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
# PLOT SETTINGS
# ============================================================

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

# ============================================================
# USER SETTINGS
# ============================================================

wss = 0                  # wall shear stress: 0, 2, 5, or 20
realization = '40'        # realization to plot

delta_t = 0.1            # s
t_average = 15 * 60      # averaging starts at 15 min

Np = 300
Lx = 0.12
Ly = 1e-3
max_area = Lx * Ly
max_len = 4.04e-5
alpha = None

# ============================================================
# READ DATA
# ============================================================

df_pos = smut.read_pos(f'results{wss}.txt')
df_bond = smut.read_bonds(f'bonds{wss}.txt')

uf_areas, uf_clusters, uf_counts = smut.cluster_data_union_find(
    df_bond, df_pos, Np, max_area, max_len, alpha
)

# ============================================================
# EXTRACT AREA FOR SELECTED REALIZATION
# ============================================================

times = []
areas = []

for key, area in uf_areas.items():

    time_step, ri = key.split('-')

    if ri != realization:
        continue

    time_step = int(time_step)
    time = time_step * delta_t

    if isinstance(area, float) and np.isnan(area):
        continue

    times.append(time)
    areas.append(float(area))

# Sort by time
order = np.argsort(times)
times = np.asarray(times)[order]
areas = np.asarray(areas)[order]

if len(times) == 0:
    raise RuntimeError(
        f'No data found for realization {realization} at WSS = {wss}.'
    )

# ============================================================
# AREA RATIO
# ============================================================

A0 = areas[0]
area_ratio = areas / A0

# ============================================================
# PLOT AREA RATIO
# ============================================================

fig = plt.figure(figsize=(7, 4.5))
ax = fig.add_axes([0.16, 0.16, 0.80, 0.76])

ax.plot(
    times / 60,
    area_ratio,
    color='darkred',
    linewidth=1.5,
    label=f'Realization {realization}'
)

ax.axvline(
    t_average / 60,
    linestyle='--',
    color='black',
    linewidth=1.5,
    label='Averaging starts (15 min)'
)

ax.axvspan(
    t_average / 60,
    times[-1] / 60,
    color='darkred',
    alpha=0.08
)

ax.set_xlabel('Time (min)')
ax.set_ylabel(r'Area ratio ($A(t)/A_0$)')
ax.set_title(
    f'WSS = {wss} dyn/cm$^2$, realization {realization}'
)

ax.set_xlim(0, 60)
ax.set_ylim(0, 1.5)

ax.legend()
ax.grid(alpha=0.25)

fig.savefig(
    f'area_ratio_wss{wss}_realization{realization}.svg',
    format='svg',
    dpi=300
)

plt.close(fig)

# ============================================================
# EXTRACT NUMBER OF CLUSTERS FOR SELECTED REALIZATION
# ============================================================

times_n = []
numbers = []

for key, value in uf_counts.items():

    time_step, ri = key.split('-')

    if ri != realization:
        continue

    time_step = int(time_step)

    times_n.append(time_step * delta_t / 60)
    numbers.append(value[0])

# Sort by time
order = np.argsort(times_n)
times_n = np.asarray(times_n)[order]
numbers = np.asarray(numbers)[order]

if len(times_n) == 0:
    raise RuntimeError(
        f'No cluster-number data found for realization '
        f'{realization} at WSS = {wss}.'
    )

# ============================================================
# PLOT NUMBER OF CLUSTERS
# ============================================================

#fig, ax = plt.subplots(figsize=(7, 4.5))
fig = plt.figure(figsize=(7, 4.5))
ax = fig.add_axes([0.16, 0.16, 0.80, 0.76])


ax.plot(
    times_n,
    numbers,
    color='navy',
    linewidth=1.5,
    label=f'Realization {realization}'
)

# 15-min averaging window
ax.axvline(
    t_average / 60,
    color='black',
    linestyle='--',
    linewidth=1.5,
    label='Averaging starts (15 min)'
)

ax.axvspan(
    t_average / 60,
    times_n[-1],
    color='navy',
    alpha=0.08
)

ax.set_xlabel('Time (min)')
ax.set_ylabel('Number of clusters')
ax.set_title(
    f'WSS = {wss} dyn/cm$^2$, realization {realization}'
)

ax.set_xlim(0, 60)
ax.set_ylim(0, 80)

ax.legend()
ax.grid(alpha=0.25)

#plt.tight_layout()


fig.savefig(
    f'cluster_number_wss{wss}_realization{realization}.svg',
    format='svg',
    dpi=300
)
plt.close(fig)


