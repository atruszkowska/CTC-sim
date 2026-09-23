import numpy as np
from scipy.spatial import ConvexHull, QhullError
import matplotlib.pyplot as plt 

def get_interparticle_distance(points):
    """ 
    Compute the interparticle distance between two particles
    if only two particles present
    """
    return np.sqrt((points[0,0]-points[1,0])*(points[0,0]-points[1,0]) + (points[0,1]-points[1,1])*(points[0,1]-points[1,1]))


def get_area(points, N_prt):

    """

    Compute the area spanned by cells using convex hull - if the cells
    are aligned, perturb their positions slightly to allow the
    computation - this does not change the actual positions
    only the input to convex hull is changed

    """
    
    # Compute the area spanned by points using convex hull
    while True:
        try:
            hull = ConvexHull(points)
            return hull.volume
    
        except:
            # Perturb the positions if all points have the same x and y coordinates
            points[:, 0] *= (1 + np.random.uniform(0, 1, size=points.shape[0])/100)
            points[:, 1] *= (1 + np.random.uniform(0, 1, size=points.shape[0])/100)


def get_fold_change_in_area(data, N_prt, N_times, N_rep):

    area_change = np.zeros((N_rep, N_times))
    cluster_count = 0

    for ir in range(N_rep):

        # Dataset portion for this repetition
        idx0 = ir * N_prt * N_times
        stride = N_prt

        # Compute the initial area
        points_initial = data[idx0:idx0 + stride, :2]
        if N_prt > 2:
            A0 = get_area(points_initial, N_prt) 
        else:
            A0 = get_interparticle_distance(points_initial)
        area_change[ir, 0] = 1

        # Compute all remaining areas, store each ratio
        for jt in range(1, N_times):
            idx = idx0 + jt * N_prt
            points = data[idx:idx + stride, :2]
            if N_prt > 2: 
                A = get_area(points, N_prt)
            else: 
                A = get_interparticle_distance(points)
            area_change[ir, jt] = A / A0
        
        if area_change[ir, -1] <= 0.3:
            cluster_count += 1

    return area_change, cluster_count

# Dataset
particle_positions = np.loadtxt('results.txt')

area_change, cluster_count = get_fold_change_in_area(data=particle_positions, N_prt=5, N_times=5, N_rep=1000)

# Save results to a file if needed
np.savetxt('area_change.txt', area_change)



