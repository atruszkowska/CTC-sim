# # # # # # # # # # # # # # # # # #
#
# Module with visualization tools 
#
# # # # # # # # # # # # # # # # # #

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle, Circle
from matplotlib import patheffects as pe
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from collections import defaultdict
import math
import matplotlib as mpl
from collections import Counter

def spy_plot(data, Nr, Ncol):
    ''' Plot sparsity pattern of 
            each sublist in data in a NrxNcol
            subplot composition; plots the data
            in columns, then progresses to rows i.e.
            Row 1: data 1, data 2, data 3;
            Row 2: data 4, data 5, data 6 '''
    
    fig, axs = plt.subplots(Nr, Ncol)
    if Nr == 1 and Ncol > 1:
        # This is separated so there is no indexing
        # error when accessing axs
        for ind, subset in enumerate(data):
            ax = axs[ind]
            ax.spy(subset, markersize=5)
    elif Nr == 1 and Ncol == 1:
        # Just one plot
        axs.spy(data[0], markersize=5)
    else:
        ind = 0
        for ir in range(0, Nr):
            for ic in range(0, Ncol):
                ax = ax[ir,ic]
                ax.spy(data[ind], markersoze=5)
                ind += 1

    plt.show()

def plot_2D_with_an(x, y_num, y_an):
    ''' Plot of numerical and analytical solutions '''

    plt.plot(x, y_num, 'b', label = 'Numerical solution')
    plt.plot(x, y_an, 'ro', label = 'Analytical solution')
    plt.xlabel('Channel width', fontsize = 18)
    plt.ylabel('Velocity', fontsize = 18)
    plt.legend(loc = 2)
    plt.grid()
    plt.show()

def plot_3D_flat(nrow, ncol, sol):
    ''' Plot a 3D surface plot, initially shown as flat 2D surface '''  

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    # Grid
    X = np.arange(0, nrow)
    Y = np.arange(0, ncol)
    X, Y = np.meshgrid(X, Y)

    # Plot
    surf = ax.plot_surface(X, Y, sol, cmap=cm.jet,
                       linewidth=0, antialiased=False)

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    ax.view_init(90, 90)
    plt.show()

def plot_velocity_field(nrow, ncol, ux, uy, geom):
    ''' Plot a velocity field as streamlines '''

    # Grid
    X = np.arange(0, nrow)
    Y = np.arange(0, ncol)
    X, Y = np.meshgrid(X, Y)

    # Velocity field components
    U = ux
    V = uy 

    # Mask - solid objects
    mask = np.logical_not(geom)
    U = np.ma.array(U, mask=mask)
    V = np.ma.array(V, mask=mask)

    # Varying color along a streamline
    fig = plt.figure(figsize=(5,5))
    ax = plt.subplot()
    strm = ax.streamplot(X, Y, U, V, linewidth=2, color='b')

    ax.imshow(~mask, alpha=0.5, cmap='gray', aspect='auto')
    ax.set_aspect('equal')
    ax.invert_yaxis()

    plt.show()

def visualize_new_heading(filename):
    headings = pd.read_csv(filename, header=None)

    headings_T = headings.T

    plt.plot(headings_T)
    plt.xlabel("Index")
    plt.ylabel("Values")
    plt.title("Visualization of new particle headings")
    plt.show()

def compare_particle_positions(file1, file2):
    df_0 = pd.read_csv(file1)
    df_1 = pd.read_csv(file2)
    
    df_0.columns = df_0.columns.str.strip()
    df_1.columns = df_1.columns.str.strip()

    pos_x_0, pos_y_0, head_0, pid_0 = df_0['position_x'], df_0['position_y'], df_0['heading'], df_0['particle_id']
    pos_x_1, pos_y_1, head_1, pid_1 = df_1['position_x'], df_1['position_y'], df_1['heading'], df_1['particle_id']
    
    u_0 = np.cos(np.radians(head_0))
    v_0 = np.sin(np.radians(head_0))

    u_1 = np.cos(np.radians(head_1))
    v_1 = np.sin(np.radians(head_1))
    
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))
    
    for particle_id in pid_0.unique():
        mask = pid_0 == particle_id
        color = plt.cm.viridis(particle_id / len(pid_0))
        axes[0].quiver(pos_x_0[mask], pos_y_0[mask], u_0[mask], v_0[mask], color = color, label=f'Particle {particle_id}', alpha=0.7)
    
    axes[0].set_xlabel('Position x')
    axes[0].set_ylabel('Position y')
    axes[0].set_title('Zero particle propagation vectors')
    
    for particle_id in pid_1.unique():
        mask = pid_1 == particle_id
        color = plt.cm.magma(particle_id / len(pid_1))
        axes[1].quiver(pos_x_1[mask], pos_y_1[mask], u_1[mask], v_1[mask], color=color, label=f'Particle {particle_id}', alpha=0.7)

    axes[1].set_xlabel('Position x')
    axes[1].set_ylabel('Position y')
    axes[1].set_title('Nonzero particle propagation vectors')

    #plt.legend()

    plt.show()

def plot_clusters(pos_df, clusters_dict, it, ir):

    key = f"{it}-{ir}"
    clusters = clusters_dict.get(key)

    if clusters is None or not clusters or not clusters[0]:
        print(f"No cluster data found for it={it}, ir={ir}")
        return

    if 'pid' not in pos_df.columns:
        pos_df = pos_df.copy()
        pos_df['pid'] = pos_df.groupby(['it', 'ir']).cumcount()
    

    # Extract positions at time t, rep ir
    p0 = pos_df[(pos_df['it'] == it) & (pos_df['ir'] == ir)].set_index("pid")

    fig, ax = plt.subplots(figsize=(8,5))
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(title or f"Clusters at it={it}, ir={ir} — {len(clusters[0])} clusters")

    colors = cm.get_cmap("tab10", 10)

    for i, comp in enumerate(clusters[0]):
        coords = p0.loc[comp, ['x', 'y']].to_numpy()
        ax.scatter(coords[:, 0], coords[:, 1],
                   s=100, color=colors(i % 10),
                   label=f"Cluster {i+1}", alpha=0.6, edgecolors='k')

    ax.legend(loc='best', fontsize='small', frameon=False)
    ax.grid(True)
    plt.tight_layout()
    plt.show()

def plot_lj_double(fname, labels):

    data = np.loadtxt(fname, skiprows=1)
    
    H = 2e-4
    sigma = 10.1
    y = data[:, 0] * 1e6

    mask = y > 0.01 * sigma
    y = y[mask]

    for k, label in enumerate(labels, 1):
        plt.plot(y, data[mask, k], lw=1.8, label=label)
    
   
    micron = "\u00B5"
    epsilon = "\u03B5"

    plt.axvline(0, color='gray', ls='--', lw=2)
    plt.axvline(H*1e6, color='gray', ls='--', lw=2)
    plt.axhline(0, color='k', ls='--', lw=0.8)
    plt.xlabel(f"Channel width, {micron}m", fontsize=20)
    plt.ylabel("Potential, J", fontsize=20)
    plt.xlim(-5, (H*1e6)+5)
    plt.ylim(-1.5, 6)

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

#    plt.title(fname)
    plt.legend()
    plt.legend(fontsize=20)

    plt.grid(True, ls=':', lw=0.5)
    plt.tight_layout()
    plt.show()


def plot_lj_single(fname, labels):

    data = np.loadtxt(fname, skiprows=1)
    
    sigma = 10.1
    y_max = 3 * sigma
    y = data[:, 0] * 1e6
    
    mask = y > 0.01 * sigma
    y = y[mask]

    for k, label in enumerate(labels, 1):
        plt.plot(y, data[mask, k], lw=1.8, label=label)

    micron = "\u00B5"
    epsilon = "\u03B5"

    plt.axvline(0, color='gray', ls='--', lw=1.2)
    plt.axhline(0, color='k', ls='--', lw=0.8)
    plt.xlabel(f"Gap y ({micron})")
    plt.ylabel("Potential V (units of {epsilon})")
    plt.title(fname)
    plt.ylim(-1.8, 6)
    #plt.xlim(-0.5, y_max*1e6)
    plt.legend()
    plt.grid(True, ls=':', lw=0.5)
    plt.tight_layout()
    plt.show()

def plot_lj_offset(fname, deltas):

    data = np.loadtxt(fname, skiprows=1)

    sigma = 10.1
    y_max = 3 * sigma
    y = data[:, 0] * 1e6

    mask = y > 0.01 * sigma
    y = y[mask]

    for idx, delta in enumerate(deltas):
        plt.plot(y, data[mask, idx+1], lw=1.8, label=fr"$\Delta = {delta:g}\, \sigma$")

    micron = "\u00B5"
    epsilon = "\u03B5"

    plt.axvline(0, color='gray', ls='--', lw=1.2)
    plt.axhline(0, color='k', ls='--', lw=0.8)
    plt.xlabel(f"Gap y ({micron})")
    plt.ylabel("Potential V (units of {epsilon})")
    plt.title(fname)
    plt.ylim(-1.8, 6)
    #plt.xlim(-0.5, y_max*1e6)
    plt.legend()
    plt.grid(True, ls=':', lw=0.5)
    plt.tight_layout()
    plt.show()

def plot_lj_eps(fname, eps_list):

    data = np.loadtxt(fname, skiprows=1)

    sigma = 10.1
    y_max = 3 * sigma
    y = data[:, 0] * 1e6

    mask = y > 0.01 * sigma
    y = y[mask]
    epsilon = "\u03B5"
    micron = "\u00B5"
    for k, label in enumerate(eps_list, 1):
        plt.plot(y, data[mask, k], lw=1.8, label=f"{epsilon} = {label}")

    plt.axvline(0, color='gray', ls='--', lw=1.2)
    plt.axhline(0, color='k', ls='--', lw=0.8)
    plt.xlabel(f"Gap y ({micron})")
    plt.ylabel("Potential V (units of {epsilon})")
    plt.title(fname)
    max_eps = max(eps_list)
    plt.ylim(-2.0*max_eps, 6*max_eps)
    #plt.xlim(-0.5, y_max*1e6)
    plt.legend()
    plt.grid(True, ls=':', lw=0.5)
    plt.tight_layout()
    plt.show()


def plot_lj_sigma(fname, sigmalist):


    data = np.loadtxt(fname, skiprows=1)

    y = data[:, 0] * 1e6
    y_max = 3*max(sigmalist)
    sig = "\u03C3"
    micron = "\u00B5"
    epsilon = "\u03B5"
    
    for idx, sigma in enumerate(sigmalist, 1):
        mask = y > 0.01 * sigma * 1e6    
        plt.plot(y[mask], data[mask, idx], lw=1.8, label=f"{sig} = {sigma*1e6:.1f} {micron}m")


    plt.axvline(0, color='gray', ls='--', lw=1.2)
    plt.axhline(0, color='k', ls='--', lw=0.8)
    plt.xlabel(f"Gap y ({micron})")
    plt.ylabel("Potential V (units of {epsilon})")
    plt.title(fname)
    plt.ylim(-1.2, 6)
    plt.xlim(-0.5, y_max*1e6)
    plt.legend()
    plt.grid(True, ls=':', lw=0.5)
    plt.tight_layout()
    plt.show()


def pick_repir(uf_areas, tol=0.01):
    """
    uf_areas: dict "it-ir" -> area
    tol: absolute tolerance on Af/A0 around the mean 
    Returns: rep_ir, mean_ratio, per_run_ratio (dict ir -> Af/A0)
    """
    # Group by run, record areas by time
    by_run = defaultdict(dict)
    for key, A in uf_areas.items():
        it_str, ir_str = key.split('-', 1)
        try:
            t = int(it_str) if it_str.isdigit() else float(it_str)
        except:
            t = it_str
        by_run[str(ir_str)][t] = A

    # Compute final ratio Af/A0 per run
    per_run_ratio = {}
    for ir, tmap in by_run.items():
        if not tmap:
            continue
        t0 = min(tmap.keys()); tf = max(tmap.keys())
        A0, Af = tmap.get(t0), tmap.get(tf)
        if any(v is None or (isinstance(v, float) and math.isnan(v)) for v in (A0, Af)) or A0 == 0:
            continue
        per_run_ratio[ir] = float(Af)/float(A0)

    if not per_run_ratio:
        raise ValueError("No usable runs to compute Af/A0.")

    mean_ratio = float(np.mean(list(per_run_ratio.values())))
    lo, hi = mean_ratio - tol, mean_ratio + tol

    # Candidates within tolerance window; if none, pick closest
    in_band = [ir for ir, r in per_run_ratio.items() if lo <= r <= hi]
    rep_ir = in_band[0] if in_band else min(per_run_ratio, key=lambda ir: abs(per_run_ratio[ir]-mean_ratio))
    return rep_ir, mean_ratio, per_run_ratio

def _components(uf_clusters, t, ir):
    comps = uf_clusters.get(f"{t}-{ir}", [])
    return comps[0] if comps and isinstance(comps[0], (list, tuple)) else comps

def _pid2comp(comps):
    m = {}
    for k, comp in enumerate(comps):
        for pid in comp:
            m[int(pid)] = k
    return m

def _translate_components_to_frame(comps_t, ids_sorted):
    """
    Ensure comps_t pids match ids_sorted (the pids we plot).
    Heuristic:
      - If max pid in comps_t < len(ids_sorted), we assume comps_t uses
        *local* 0..n-1 indices, so we translate via inv_local.
      - Otherwise we assume comps_t already uses the same ids as ids_sorted.
    """
    if not comps_t:
        return comps_t

    max_pid = max(max(comp) for comp in comps_t if comp)
    n = len(ids_sorted)
    if max_pid < n:
        # comps use local indices 0..n-1 -> translate to ids_sorted values
        inv_local = {i: pid for i, pid in enumerate(ids_sorted)}
        comps_tx = [[inv_local[i] for i in comp if i in inv_local] for comp in comps_t]
    else:
        # assume they already use the same pid space as ids_sorted
        comps_tx = comps_t
    return comps_tx


def _labels_for_time(ids_sorted, comps_t, base_pid2comp=None, next_color=0):
    """
    Make per-particle labels aligned with ids_sorted, assigning the same label
    to *every* member of a cluster. Singletons get -1.
    """
    pid2lab = {}

    for comp in comps_t:
        inherited = None
        if base_pid2comp:
            votes = [base_pid2comp.get(pid) for pid in comp if pid in base_pid2comp]
            if votes:
                counts = {}
                for v in votes:
                    counts[v] = counts.get(v, 0) + 1
                inherited = max(counts, key=counts.get)

        lab = inherited if inherited is not None else next_color
        if inherited is None:
            next_color += 1

        for pid in comp:
            pid2lab[pid] = lab

    labels = [pid2lab.get(pid, -1) for pid in ids_sorted]
    return labels, next_color

def plot_cells_clusters(df_pos, uf_clusters, rep_ir, times_to_plot, Lx, Ly,
                       particle_radius, center_channel=True,
                       radius_mode="screen", 
                       size_pts=10,          
                       vis_scale=2.0):       
    """
    Visualize particle clusters in a long microfluidic channel at selected times.

    The function joins cluster membership (from `uf_clusters`) to particle
    positions (from `df_pos`) by pid, assigns a cluster label to every member,
    and plots each particle as a disk. Particles sharing the same label are
    rendered with the same color; singletons (label = -1) are drawn in grey.

    Parameters
    ----------
    df_pos : pandas.DataFrame
        Particle positions. Must include columns:
        - 'x', 'y' : coordinates in channel units
        - 'it'     : integer time index
        - 'ir'     : run index (string/int)
        - 'id'     : stable particle id (pid) consistent with `uf_clusters`
          (if your pid column is named differently, rename it to 'id' first).

    uf_clusters : dict
        Mapping from (time, run) → list of components.
        Each component is a list of pids. If components are stored as local
        0..n-1 indices per frame, they are translated to the frame’s pids
        internally before labeling.

    rep_ir : str or int
        Run index to visualize (matched against df_pos['ir']).

    times_to_plot : list[int]
        Time indices to render, in order (e.g., [0, 6, 10]).

    Lx, Ly : float
        Channel length and height (same units as 'x'/'y').

    particle_radius : float
        Physical particle radius in data units (same units as 'y').

    center_channel : bool, default True
        If True, plot y in [-Ly/2, +Ly/2]; otherwise in [0, Ly].

    radius_mode : {"screen", "data"}, default "screen"
        Controls how marker size is specified:
        - "screen": size is constant in screen points (use `size_pts`).
        - "data"  : size uses physical data units (uses `particle_radius`
          multiplied by `vis_scale` for visibility).

    size_pts : float, default 10
        Marker size in points (used only when `radius_mode=="screen"`).
        Matplotlib's scatter uses area; we pass `size_pts**2`.

    vis_scale : float, default 2.0
        Visual multiplier for `particle_radius` when `radius_mode=="data"`.
        Use a modest factor (e.g., 1.5–3) so overlapping particles remain
        distinguishable without inflating clusters.

    Notes
    -----
    - Colors are assigned per-cluster label with a cyclic tab20 colormap.
    - Axis limits are fixed to the domain: x ∈ [0, Lx], y ∈ [y_lo, y_hi].
    - A rectangle outlines the channel; an image 'cluster_plot.png' is saved.

    Returns
    -------
    None
        Displays the figure and saves 'cluster_plot.png'.
    """
    df_pos = df_pos.sort_values(['ir', 'it']).reset_index(drop=True)
    
    if 'id' not in df_pos.columns:
        df_pos['id'] = df_pos.groupby(['it', 'ir']).cumcount().astype(int)
    df_pos['it'] = pd.to_numeric(df_pos['it'], errors='coerce').astype('Int64')
    df_pos['ir'] = df_pos['ir'].astype(str)

    if center_channel:
        df_pos = df_pos.copy()
        df_pos['y_plot'] = df_pos['y'] - Ly/2.0
        ycol = 'y_plot'; y_lo, y_hi = -Ly/2.0, +Ly/2.0
    else:
        ycol = 'y';      y_lo, y_hi = 0.0, Ly

    g = df_pos[df_pos['ir'] == str(rep_ir)].copy()
    cmap = plt.cm.get_cmap("tab20", 20)

    fig, axes = plt.subplots(len(times_to_plot), 1, figsize=(18, 6))
    axes = np.atleast_1d(axes)

    t0 = times_to_plot[0]
    base_comps = _components(uf_clusters, t0, rep_ir)
    base_pid2comp = _pid2comp(base_comps)
    next_color = (len(base_comps) or 0)

    for ax, t in zip(axes, times_to_plot):
        gp = g[g['it'] == t].copy()
        if gp.empty:
            ax.set_axis_off(); ax.set_title(f"t = {t} (no data)"); continue

        gp = gp.sort_values('id')
        ids_sorted = gp['id'].astype(int).tolist()
        xs = gp['x'].to_numpy(); ys = gp[ycol].to_numpy()

        comps_t_raw = _components(uf_clusters, t, rep_ir)
        comps_t = _translate_components_to_frame(comps_t_raw, ids_sorted)
        labels, next_color = _labels_for_time(ids_sorted, comps_t, base_pid2comp, next_color)
        labels = np.asarray(labels)
        print(f"t={t}: cluster sizes (by labels) ->", Counter([l for l in labels if l >= 0]))
        # frame
        ax.add_patch(Rectangle((0, y_lo), Lx, y_hi-y_lo, fill=False,
                               linewidth=1.0, edgecolor=(0,0,0,0.7), zorder=1))

        # draw singletons (grey, under)
        sing = labels < 0
        if radius_mode == "screen":
            ax.scatter(xs[sing], ys[sing], s=size_pts**2, c=[(0.5,0.5,0.5,0.3)],
                       edgecolors=(0,0,0,0.25), linewidths=0.4, zorder=2)
        else:
            r_plot = vis_scale * float(particle_radius)
            for x, y in zip(xs[sing], ys[sing]):
                ax.add_patch(Circle((x, y), r_plot,
                                    facecolor=(0.5,0.5,0.5,0.3),
                                    edgecolor=(0,0,0,0.25), linewidth=0.4, zorder=2))

        # clustered cells (colored, thin white edge so overlaps are visible)
        cl = np.where(labels >= 0)[0]
        if radius_mode == "screen":
            colors = [cmap(int(labels[i]) % 20) for i in cl]
            ax.scatter(xs[cl], ys[cl], s=(size_pts**2),
                       c=colors, edgecolors='white', linewidths=0.6, zorder=3)
        else:
            r_plot = vis_scale * float(particle_radius)
            for i in cl:
                face = cmap(int(labels[i]) % 20)
                circ = Circle((xs[i], ys[i]), r_plot,
                              facecolor=face, edgecolor='white',
                              linewidth=0.6, zorder=3)
                ax.add_patch(circ)

        # domain-aligned limits
        ax.set_xlim(0.0, Lx)
        ax.set_ylim(y_lo, y_hi)
        ax.set_aspect('auto')     
        ax.set_axis_off()
        ax.set_title(f"t = {t}s")

    plt.tight_layout()
    plt.savefig("cluster_plot.png", dpi=300, bbox_inches="tight")
    plt.show()

