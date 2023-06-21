import numpy as np
from numba import njit


# Equivalent to np.linalg.norm([a, b], axis=0), which Numba has not implemented
@njit
def _norm(a, b):
    return np.sqrt((a**2) + (b**2))


# Calculate horizontal distance vector
@njit
def _horizontal_distance(x, y):
    return np.hstack((np.array([0]), np.cumsum(_norm(np.diff(x), np.diff(y)))))


# Compute along-path length
@njit
def _path_length(d, z):
    return _norm(np.diff(d), np.diff(z)).sum()


# Compute DIRECT path
@njit
def _direct_path(d, z):
    return (z[-1] - z[0]) / (d[-1] - d[0]) * (d - d[0]) + z[0]  # d - d[0] for intercept


# [Recursively] compute SHORTEST DIFFRACTED path
@njit
def _shortest_diffracted_path(d, z):
    # Compute the direct path
    direct_path = _direct_path(d, z)

    # If z is everywhere "on" or "underneath" the direct path, we're done (first
    # we mask values that are ~equal; then we check for "less than")
    isclose = np.abs(z - direct_path) < np.finfo(np.float32).eps
    if (z[~isclose] < direct_path[~isclose]).all():
        return direct_path

    # Location of maximum of profile (detrended using direct_path)
    max_ind = np.argmax(z - direct_path)

    # Split profile into two pieces here (including common midpoint in both)
    left = slice(None, max_ind + 1)
    right = slice(max_ind, None)

    # Recursively call this function
    path_left = _shortest_diffracted_path(d[left], z[left])
    path_right = _shortest_diffracted_path(d[right], z[right])

    # Join at common midpoint, removing duplicate
    return np.hstack((path_left[:-1], path_right))
