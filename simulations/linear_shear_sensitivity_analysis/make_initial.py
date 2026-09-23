import numpy as np
import math


# ============================================================
# USER PARAMETERS
# ============================================================

# Total number of cells
n_cells = 300

# Number of clusters of each size
n_single = int(56)#*2/3) + 9
n_2cluster =int( 25)#*2/3)
n_3cluster =int( 16)#*2/3)
n_4cluster =int( 12)#*2/3)
n_5cluster =int( 10)#*2/3)
n_6cluster =int( 8)#*2/3)

# Cell radius (same for every cell)
cell_radius = 1.01e-5

# Maximum center-to-center distance for bonding
max_bond_distance = 2.25*cell_radius

# Rectangular initialization domain
xmin = 0.0
xmax = 0.12
ymin = 0.0
ymax = 1e-3

# Output file
output_file = "initial_cells.txt"

# Maximum number of attempts to place each cluster
max_attempts = 100000

# ============================================================
# CHECK INPUT
# ============================================================

calculated_cells = (
    n_single
    + 2 * n_2cluster
    + 3 * n_3cluster
    + 4 * n_4cluster
    + 5 * n_5cluster
    + 6 * n_6cluster
)

if calculated_cells != n_cells:
    raise ValueError(
        f"\nERROR: Specified n_cells = {n_cells}, "
        f"but the requested cluster composition contains "
        f"{calculated_cells} cells."
    )

min_separation = 2.0 * cell_radius

if max_bond_distance <= min_separation:
    raise ValueError(
        "\nERROR: max_bond_distance must be greater than "
        "2 * cell_radius."
    )

if xmin >= xmax or ymin >= ymax:
    raise ValueError("\nERROR: Invalid rectangular domain.")


# ============================================================
# CLUSTER GEOMETRIES
# ============================================================

def generate_cluster_geometry(n):
    """
    Generate a compact connected geometry.

    Distances are scaled later so that neighboring cells are
    separated by more than 2R and less than max_bond_distance.
    """

    if n == 1:

        points = np.array([
            [0.0, 0.0]
        ])

    elif n == 2:

        points = np.array([
            [-0.5, 0.0],
            [ 0.5, 0.0]
        ])

    elif n == 3:

        # Triangle
        points = np.array([
            [-0.5, 0.0],
            [ 0.5, 0.0],
            [ 0.0, 0.866025]
        ])

    elif n == 4:

        # Square
        points = np.array([
            [-0.5, -0.5],
            [ 0.5, -0.5],
            [ 0.5,  0.5],
            [-0.5,  0.5]
        ])

    elif n == 5:

        # Compact cross
        points = np.array([
            [ 0.0,  0.0],
            [-1.0,  0.0],
            [ 1.0,  0.0],
            [ 0.0, -1.0],
            [ 0.0,  1.0]
        ])

    elif n == 6:

        # Two rows of three cells
        points = np.array([
            [-1.0, -0.5],
            [ 0.0, -0.5],
            [ 1.0, -0.5],
            [-1.0,  0.5],
            [ 0.0,  0.5],
            [ 1.0,  0.5]
        ])

    else:
        raise ValueError("Cluster size must be between 1 and 6.")

    return points


# ============================================================
# RANDOMIZE CLUSTER
# ============================================================

def randomize_cluster(n):

    points = generate_cluster_geometry(n)

    # --------------------------------------------------------
    # Find the smallest distance in the base geometry.
    # This corresponds to the nearest-neighbor spacing.
    # --------------------------------------------------------

    min_dist = np.inf

    for i in range(n):
        for j in range(i + 1, n):

            distance = np.linalg.norm(
                points[i] - points[j]
            )

            min_dist = min(min_dist, distance)

    # Choose a target bond distance safely between the
    # non-overlap limit and the maximum bonding distance.
    target_distance = (
        min_separation + max_bond_distance
    ) / 2.0

    scale = target_distance / min_dist

    points *= scale

    # --------------------------------------------------------
    # Random rotation
    # --------------------------------------------------------

    theta = np.random.uniform(
        0.0,
        2.0 * math.pi
    )

    rotation = np.array([
        [math.cos(theta), -math.sin(theta)],
        [math.sin(theta),  math.cos(theta)]
    ])

    points = points @ rotation.T

    # --------------------------------------------------------
    # Small random perturbation
    # --------------------------------------------------------

    perturbation_size = 0.10 * (
        max_bond_distance - min_separation
    )

    perturbations = np.random.uniform(
        -perturbation_size,
        perturbation_size,
        size=(n, 2)
    )

    points += perturbations

    return points


# ============================================================
# CHECK CLUSTER CONNECTIVITY
# ============================================================

def cluster_is_connected(cluster):

    n = len(cluster)

    # Build the bond graph
    connected = [False] * n

    connected[0] = True

    changed = True

    while changed:

        changed = False

        for i in range(n):

            if not connected[i]:
                continue

            for j in range(n):

                if connected[j]:
                    continue

                distance = np.linalg.norm(
                    cluster[i] - cluster[j]
                )

                if (
                    distance > min_separation
                    and distance < max_bond_distance
                ):
                    connected[j] = True
                    changed = True

    return all(connected)


# ============================================================
# CHECK CLUSTER VALIDITY
# ============================================================

def cluster_is_valid(cluster, existing_cells):

    n = len(cluster)

    # --------------------------------------------------------
    # 1. No overlap within cluster
    # --------------------------------------------------------

    for i in range(n):

        for j in range(i + 1, n):

            distance = np.linalg.norm(
                cluster[i] - cluster[j]
            )

            if distance <= min_separation:
                return False

    # --------------------------------------------------------
    # 2. Cluster must be connected by bonds
    # --------------------------------------------------------

    if n > 1:

        if not cluster_is_connected(cluster):
            return False

    # --------------------------------------------------------
    # 3. No overlap or unintended bond with existing cells
    # --------------------------------------------------------

    for new_cell in cluster:

        for old_cell in existing_cells:

            distance = np.linalg.norm(
                new_cell - old_cell
            )

            # No overlap
            if distance <= min_separation:
                return False

            # No bond between separate clusters
            if distance < max_bond_distance:
                return False

    return True


# ============================================================
# PLACE CLUSTER
# ============================================================

def place_cluster(n, existing_cells):

    for attempt in range(max_attempts):

        cluster = randomize_cluster(n)

        # ----------------------------------------------------
        # Determine allowable center location
        # ----------------------------------------------------

        cluster_x_min = np.min(cluster[:, 0])
        cluster_x_max = np.max(cluster[:, 0])
        cluster_y_min = np.min(cluster[:, 1])
        cluster_y_max = np.max(cluster[:, 1])

        center_x_min = xmin - cluster_x_min
        center_x_max = xmax - cluster_x_max

        center_y_min = ymin - cluster_y_min
        center_y_max = ymax - cluster_y_max

        if center_x_min >= center_x_max:
            continue

        if center_y_min >= center_y_max:
            continue

        # ----------------------------------------------------
        # Random cluster center
        # ----------------------------------------------------

        center_x = np.random.uniform(
            center_x_min,
            center_x_max
        )

        center_y = np.random.uniform(
            center_y_min,
            center_y_max
        )

        cluster[:, 0] += center_x
        cluster[:, 1] += center_y

        # ----------------------------------------------------
        # Check validity
        # ----------------------------------------------------

        if cluster_is_valid(
            cluster,
            existing_cells
        ):
            return cluster

    raise RuntimeError(
        f"\nERROR: Could not place a {n}-cell cluster "
        f"after {max_attempts} attempts.\n"
        f"Try increasing the domain size or reducing "
        f"the number of cells/clusters."
    )


# ============================================================
# GENERATE INITIAL CONFIGURATION
# ============================================================

cells = []

cluster_sizes = (
    [1] * n_single
    + [2] * n_2cluster
    + [3] * n_3cluster
    + [4] * n_4cluster
    + [5] * n_5cluster
    + [6] * n_6cluster
)

# Randomize placement order
np.random.shuffle(cluster_sizes)

print(cluster_sizes)


for cluster_size in cluster_sizes:

    cluster = place_cluster(
        cluster_size,
        cells
    )

    for cell in cluster:
        cells.append(cell)


# ============================================================
# FINAL CHECK
# ============================================================

if len(cells) != n_cells:

    raise RuntimeError(
        f"ERROR: Generated {len(cells)} cells instead of "
        f"{n_cells}."
    )

cells = np.array(cells)

# Check domain boundaries
if (
    np.any(cells[:, 0] < xmin)
    or np.any(cells[:, 0] > xmax)
    or np.any(cells[:, 1] < ymin)
    or np.any(cells[:, 1] > ymax)
):

    raise RuntimeError(
        "ERROR: One or more cells are outside the domain."
    )


# ============================================================
# WRITE OUTPUT
# ============================================================

np.savetxt(
    output_file,
    cells,
    fmt="%.8f",
    delimiter=" "
)


# ============================================================
# SUMMARY
# ============================================================

print(f"Generated {len(cells)} cells.")
print(f"  Singles:     {n_single}")
print(f"  2-cell:      {n_2cluster}")
print(f"  3-cell:      {n_3cluster}")
print(f"  4-cell:      {n_4cluster}")
print(f"  5-cell:      {n_5cluster}")
print(f"  6-cell:      {n_6cluster}")
print(f"\nInitial configuration written to: {output_file}")
