import pandas as pd
import numpy as np
import math
from typing import Dict, List, Tuple, Optional
from scipy.spatial import ConvexHull, QhullError
import alphashape  
from shapely.geometry import Polygon, MultiPolygon

# Read particle coordinates into a data frame 
def read_pos(pos_file: str):
    """
    Read the positions and return a DataFrame:
        - pos has columns ['x','y','it','ir']
    """
    # Read raw, whitespace-delimited, no header
    pos   = pd.read_csv(pos_file,  delim_whitespace=True, header=None,
                        comment='#', on_bad_lines='skip')

    # Assign pos columns
    if pos.shape[1] == 4:
        pos.columns = ["x", "y", "it", "ir"]
    else:
        raise ValueError(f"Unexpected # of columns in pos_file: {pos.shape[1]}")

    return pos

# Read IDs of bonded particles into a data frame 
def read_bonds(bond_file: str):
    """
    Read the bond IDs and return a DataFrame:
      - bonds has columns ['ID-1','ID-2','it','ir']
      - ID-1 is always smaller than ID-2
    """
    # Read raw, whitespace-delimited, no header
    bonds = pd.read_csv(bond_file, delim_whitespace=True, header=None,
                        comment='#', on_bad_lines='skip')

    # Assign bonds columns
    if bonds.shape[1] == 3:
        bonds.columns = ["ij", "it", "ir"]
    else:
        raise ValueError(f"Unexpected # of columns in bond_file: {bonds.shape[1]}")

    # Split the first column into two
    temp_bonds = pd.DataFrame(bonds.ij.str.split('-').tolist(), columns = ['ID-1','ID-2'])
    
    bonds = pd.concat([temp_bonds, bonds.iloc[:, 1:]], axis=1)

    return bonds

# Union find for getting the clusters from bonded particles
class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0]*n

    def find(self, x):
        # path-compression
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return
        # union by rank
        if self.rank[rx] < self.rank[ry]:
            self.parent[rx] = ry
        elif self.rank[ry] < self.rank[rx]:
            self.parent[ry] = rx
        else:
            self.parent[ry] = rx
            self.rank[rx] += 1

# Cluster compositions, numbers, and areas using union find
def cluster_data_union_find(bonds: pd.DataFrame, pos: pd.DataFrame, num_particles: int, max_area: float, max_len: float, alpha: Optional[float] = None):

    """
    bonds: DataFrame with columns ['ID-1','ID-2','it','ir']
    pos:   DataFrame with columns ['x','y','it','ir']
    num_particles: number of particles per repetition
    max_area: correction for PBC, if a cluster area is larger than this - discard it
    max_len: correction for PBC, if a distance between two particle cluster larger than this - discard it
    alpha: alpha-shape parameter

    Returns:
        areas (dict): Total area of all clusters at each collected time and all repetitions
        cluster_lists (dict): key is it-ir and value is a list of clusters
        cluster_counts (dict): key is it-ir and value is a list of #clusters
    """

    if 'pid' not in pos.columns:
        pos = pos.copy()
        pos['pid'] = pos.groupby(['it','ir']).cumcount()

    if alpha is None:
        spans = pos[['x','y']].agg(lambda col: col.max() - col.min())
        alpha = 0.1 * spans.mean()

    areas = {}
    cluster_lists = {}
    cluster_counts = {}

    for it, df_it in pos.groupby('it'):

        for ir, pos_grp in df_it.groupby('ir'):
            areas_it, lists_it, counts_it = [], [], []

            p = pos_grp.copy()
            uf = UnionFind(num_particles)

            rep_bonds = bonds[(bonds["it"] == it) & (bonds["ir"] == ir)]
            for i, j in zip(rep_bonds["ID-1"], rep_bonds["ID-2"]):
                uf.union(int(i), int(j))

            clusters = {}
            for pid in range(num_particles):
                root = uf.find(pid)
                clusters.setdefault(root, []).append(pid)

            area = 0.0
            multi = []
            for comp in clusters.values():
                if len(comp) < 2:
                    # Not a cluster
                    #areas.append(0.0)
                    continue
                else:
                    mask = p["pid"].isin(comp)
                    pts = p.loc[mask, ["x", "y"]].to_numpy()
                    cluster_area = get_cluster_area(pts, alpha, 10)
                    if (len(comp) > 2 and cluster_area <= max_area) or (len(comp) == 2 and cluster_area <= max_len):
                        area += cluster_area
                        multi.append(comp)
  
            lists_it.append(multi)
            counts_it.append(len(lists_it[-1]))
                                  
            areas[str(it) + '-' + str(ir)] = area
            cluster_lists[str(it) + '-' + str(ir)] = lists_it
            cluster_counts[str(it) + '-' + str(ir)] = counts_it

    return areas, cluster_lists, cluster_counts

# Compute area using the shapes of clustered particles
def get_cluster_area(points: np.ndarray,
                     alpha: float,
                     hull_thresh: int = 10) -> float:

    """
    Compute the area of a point-cluster, using convex hull for small clusters
    and alpha-shape for larger ones. Falls back to interparticle distance for
    collinear cases.

    Parameters
    ----------
    points : (N,2) np.ndarray
        Coordinates of the cluster.
    alpha : float
        alpha parameter controlling the detail of the alpha-shape.
    hull_thresh : int, optional
        If N <= hull_thresh, use convex hull instead of alpha-shape.

    Returns
    -------
    float
        The area of the cluster.
    """

    N = points.shape[0]

    # 0 or 1 point => should not be here
    assert N >= 2

    # Exactly two points => interparticle distance
    if N == 2:
        return get_interparticle_distance(points)
    else:
        if N <= hull_thresh:
            # Small clusters => convex hull
            try:
                hull = ConvexHull(points)
                return hull.volume
            except QhullError as e:
                # Points are collinear
                print("QhullError encountered:", e) 
                return get_interparticle_distance(points)
        else:
            # Larger clusters => alpha-shape
            try:
                shape = alphashape.alphashape(points, alpha)
                if shape.geom_type == 'Polygon':
                    return shape.area
                if shape.geom_type == 'MultiPolygon':
                    return sum(poly.area for poly in shape.geoms)
            except Exception:
                return get_interparticle_distance(points)

# Get interparticle distance between two particles 
def get_interparticle_distance(points):
    """
    For a cluster with exactly two points, returns their Euclidean distance.
    """
    return np.linalg.norm(points[0] - points[1])

