import numpy as np
import xarray as xr
from pygmt.datasets import load_earth_relief
from pyproj import CRS, Transformer

# [m] Target horizontal spacing for profile (determines # points) - does not seem to
# slow down code much if this is decreased
TARGET_SPACING = 10


# Helper function to calculate horizontal difference vector for a profile DataArray
def _horizontal_distance(profile):
    return np.hstack(
        [0, np.cumsum(np.linalg.norm([np.diff(profile.x), np.diff(profile.y)], axis=0))]
    )


# Helper function for computing along-path length
def _path_length(d, z):
    return np.linalg.norm([np.diff(d), np.diff(z)], axis=0).sum()


# Helper function for computing direct path
def _direct_path(d, z):
    return (z[-1] - z[0]) / (d[-1] - d[0]) * (d - d[0]) + z[0]  # d - d[0] for intercept


# [Recursive] helper function for computing shortest diffracted path
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


def calculate_paths(src_lat, src_lon, rec_lat, rec_lon):
    """TODO: Fill in info here

    Note:
        Coordinates are expected to be in the WGS 84 datum.

    Args:
        src_lat (int or float): Source latitude
        src_lon (int or float): Source longitude
        rec_lat (int, float, list, tuple, or numpy.ndarray): One or more receiver
            latitudes
        rec_lon (int, float, list, tuple, or numpy.ndarray): One or more receiver
            longitudes

    Returns:
        TODO: Decide on output type etc.
    """

    # Type checks
    message = 'src_lat and src_lon must both be scalars!'
    assert np.isscalar(src_lat) and np.isscalar(src_lon), message

    # Type conversion, so we can iterate
    rec_lats = np.atleast_1d(rec_lat)
    rec_lons = np.atleast_1d(rec_lon)

    # Compute region based on provided source-receiver geometry (TODO: add buffer here?)
    region = [
        np.min([src_lon, rec_lons.min()]),  # xmin
        np.max([src_lon, rec_lons.max()]),  # xmax
        np.min([src_lat, rec_lats.min()]),  # ymin
        np.max([src_lat, rec_lats.max()]),  # ymax
    ]

    # Get DEM using PyGMT (using highest resolution here)
    dem = load_earth_relief(resolution='01s', region=region, use_srtm=True)
    dem.rio.write_crs(dem.horizontal_datum, inplace=True)

    # Clean DEM before going further
    dem = dem.squeeze(drop=True).rename('elevation')
    dem.attrs = {}

    # Project DEM to UTM
    utm_crs = CRS(dem.rio.estimate_utm_crs(datum_name='WGS 84'))
    dem_utm = dem.rio.reproject(utm_crs).drop('spatial_ref')

    # Get UTM coords for source and receivers
    proj = Transformer.from_crs(utm_crs.geodetic_crs, utm_crs)
    src_x, src_y = proj.transform(src_lat, src_lon)
    rec_xs, rec_ys = proj.transform(rec_lats, rec_lons)

    # Iterate over all receivers (= source-receiver pairs), computing elevation profiles
    profiles = []
    total_its = rec_lats.size
    counter = 0
    print('Computing DEM profiles...')
    for rec_x, rec_y in zip(rec_xs, rec_ys):

        # Determine # of points in profile
        dist = np.linalg.norm([src_x - rec_x, src_y - rec_y])
        n = int(np.ceil(dist / TARGET_SPACING))

        # Make profile, clean up, and add to list
        profile = dem_utm.interp(
            x=xr.DataArray(np.linspace(src_x, rec_x, n)),
            y=xr.DataArray(np.linspace(src_y, rec_y, n)),
            method='linear',
        )
        profile = profile.assign_coords(dim_0=_horizontal_distance(profile))
        profile = profile.rename(dim_0='distance')
        profiles.append(profile)

        # Print progress
        counter += 1
        print('{:.1f}%'.format((counter / total_its) * 100), end='\r')

    print('Done')

    # Iterate over all profiles, calculating paths
    ds_list = []
    for profile in profiles:

        # Ensure numpy.ndarray type for function input
        d = profile.distance.values
        z = profile.values

        # Compute DIRECT path
        direct_path = _direct_path(d, z)
        direct_path_len = _path_length(d, direct_path)

        # Compute SHORTEST DIFFRACTED path
        diff_path = _shortest_diffracted_path(d, z)
        diff_path_len = _path_length(d, diff_path)

        # Make nice Dataset of all info
        ds = xr.Dataset(
            {
                profile.name: profile,
                'direct_path': ('distance', direct_path, dict(length=direct_path_len)),
                'diffracted_path': ('distance', diff_path, dict(length=diff_path_len)),
            },
            attrs=dict(
                path_length_difference=diff_path_len - direct_path_len, units='m'
            ),
        )
        ds_list.append(ds)

    return ds_list
