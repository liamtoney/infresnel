import numpy as np


# Calculate horizontal distance vector for a profile DataArray
def _horizontal_distance(profile):
    return np.hstack(
        [0, np.cumsum(np.linalg.norm([np.diff(profile.x), np.diff(profile.y)], axis=0))]
    )


# Compute along-path length
def _path_length(d, z):
    return np.linalg.norm([np.diff(d), np.diff(z)], axis=0).sum()


# Compute DIRECT path
def _direct_path(d, z):
    return (z[-1] - z[0]) / (d[-1] - d[0]) * (d - d[0]) + z[0]  # d - d[0] for intercept


# [Recursively] compute SHORTEST DIFFRACTED path
def _shortest_diffracted_path(d, z):

    # Compute the direct path
    direct_path = _direct_path(d, z)

    # If z is everywhere "on" or "underneath" the direct path, we're done (first
    # we mask values that are ~equal; then we check for "less than")
    isclose = np.isclose(z, direct_path)
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
    return np.hstack([path_left[:-1], path_right])
