import time
from pathlib import Path

import numpy as np
import pandas as pd
import pygmt
import xarray as xr
from pyproj import CRS, Transformer
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info
from rasterio.enums import Resampling
from scipy.interpolate import RectBivariateSpline

from .helpers import (
    _direct_path,
    _horizontal_distance,
    _path_length,
    _shortest_diffracted_path,
)


def calculate_paths(
    src_lat, src_lon, rec_lat, rec_lon, dem_file=None, full_output=False
):
    """Calculate elevation profiles, direct paths, and shortest diffracted paths.

    Paths are calculated for a given DEM (either a user-supplied `dem_file` or
    automatically downloaded 1 arc-second SRTM data) and an arbitrary number of
    source-receiver pairs. By default, the function returns only the path length
    differences. If `full_output` is set to `True`, then the complete path information
    (lengths, coordinates, etc.) and the DEM used are returned.

    Note:
        Input coordinates are expected to be in the WGS 84 datum. DEM file vertical
        units are expected to be meters.

    Args:
        src_lat (int or float): Source latitude
        src_lon (int or float): Source longitude
        rec_lat (int, float, list, tuple, or :class:`~numpy.ndarray`): One or more receiver
            latitudes
        rec_lon (int, float, list, tuple, or :class:`~numpy.ndarray`): One or more receiver
            longitudes
        dem_file (str or None): Path to DEM file (if `None`, then SRTM data are used)
        full_output (bool): Toggle outputting full profile/path information and DEM, vs.
            just the path length differences

    Returns:
        If `full_output` is `False` — a :class:`~numpy.ndarray` of path length differences [m],
        one per source-receiver pair

        If `full_output` is `True` — a tuple of the form ``(ds_list, dem)`` where
        ``ds_list``  is a list of :class:`~xarray.Dataset` objects, one per source-receiver pair,
        containing full profile and path information, and ``dem`` is a :class:`~xarray.DataArray`
        containing the UTM-projected DEM used to compute the profiles
    """

    # Type checks
    message = 'src_lat and src_lon must both be scalars!'
    assert np.isscalar(src_lat) and np.isscalar(src_lon), message

    # Type conversion, so we can iterate
    rec_lats = np.atleast_1d(rec_lat)
    rec_lons = np.atleast_1d(rec_lon)

    print('Loading and projecting DEM...')
    if dem_file is not None:
        # Load user-provided DEM, first checking if it exists
        dem_file = Path(str(dem_file)).expanduser().resolve()
        assert dem_file.is_file(), 'dem_file does not exist!'
        dem = xr.open_dataarray(dem_file)
    else:
        # Get SRTM data using PyGMT, computing region (buffered by 5% in each direction)
        # based on provided source-receiver geometry (have to manually write the CRS for
        # these files)
        xmin = np.min([src_lon, rec_lons.min()])
        xmax = np.max([src_lon, rec_lons.max()])
        ymin = np.min([src_lat, rec_lats.min()])
        ymax = np.max([src_lat, rec_lats.max()])
        x_buffer = (xmax - xmin) * 0.05
        y_buffer = (ymax - ymin) * 0.05
        region = [xmin - x_buffer, xmax + x_buffer, ymin - y_buffer, ymax + y_buffer]
        with pygmt.config(GMT_VERBOSE='e'):  # Suppress warnings
            dem = pygmt.datasets.load_earth_relief(
                resolution='01s', region=region, use_srtm=True
            )
        dem.rio.write_crs(dem.horizontal_datum, inplace=True)

    # Clean DEM before going further
    dem = dem.squeeze(drop=True).rename('elevation')

    # Project DEM to UTM, further relabeling
    utm_crs = CRS(dem.rio.estimate_utm_crs(datum_name='WGS 84'))
    dem_utm = dem.rio.reproject(utm_crs, resampling=Resampling.cubic_spline)
    units = dict(units='m')
    dem_utm.attrs = units
    for coordinate in 'x', 'y':
        dem_utm[coordinate].attrs = units
    print('Done\n')

    # Determine target spacing of interpolated profiles from DEM spacing - decreasing
    # the spacing makes things slower! TODO: Is oversampling actually needed w/ spline interpolation?
    mean_resolution = np.abs(dem_utm.rio.resolution()).mean()
    target_spacing = mean_resolution / 2  # [m] Oversample to avoid aliasing
    print(
        f'DEM spacing = {mean_resolution:.2f} m -> profile spacing = {target_spacing:.2f} m\n'
    )

    # Get UTM coords for source and receivers
    proj = Transformer.from_crs(utm_crs.geodetic_crs, utm_crs)
    src_x, src_y = proj.transform(src_lat, src_lon)
    rec_xs, rec_ys = proj.transform(rec_lats, rec_lons)

    # Fit bivariate spline to DEM (slow for very high resolution DEMs!)
    print('Fitting spline to DEM...')
    x = dem_utm.x
    y = dem_utm.y
    z = dem_utm.fillna(0)  # Can't have NaNs in z
    if not pd.Series(x).is_monotonic_increasing:
        x = x[::-1]
        z = np.fliplr(z)
    if not pd.Series(y).is_monotonic_increasing:
        y = y[::-1]
        z = np.flipud(z)
    spline = RectBivariateSpline(x=x, y=y, z=z.T)  # x and y are monotonic increasing
    print('Done\n')

    # Iterate over all receivers (= source-receiver pairs), calculating paths
    ds_list = []
    counter = 0
    print(f'Computing {rec_lats.size} paths...')
    for rec_x, rec_y in zip(rec_xs, rec_ys):

        # Determine # of points in profile
        dist = np.linalg.norm([src_x - rec_x, src_y - rec_y])
        n = int(np.ceil(dist / target_spacing))

        # Make profile by evaluating spline
        xvec = np.linspace(src_x, rec_x, n)
        yvec = np.linspace(src_y, rec_y, n)
        profile = xr.DataArray(
            spline.ev(xvec, yvec),
            dims='distance',
            coords=dict(x=('distance', xvec, units), y=('distance', yvec, units)),
            attrs=units,
        )
        profile = profile.assign_coords(
            distance=_horizontal_distance(profile.x.values, profile.y.values)
        )
        profile.distance.attrs = units

        # Compute DIRECT path
        direct_path = _direct_path(profile.distance.values, profile.values)
        direct_path_len = _path_length(profile.distance.values, direct_path)

        # Compute SHORTEST DIFFRACTED path
        diff_path = _shortest_diffracted_path(profile.distance.values, profile.values)
        diff_path_len = _path_length(profile.distance.values, diff_path)

        # Make nice Dataset of all info
        ds = xr.Dataset(
            {
                'elevation': profile,
                'direct_path': (
                    'distance',
                    direct_path,
                    dict(length=direct_path_len, **units),
                ),
                'diffracted_path': (
                    'distance',
                    diff_path,
                    dict(length=diff_path_len, **units),
                ),
            },
            attrs=dict(path_length_difference=diff_path_len - direct_path_len, **units),
        )
        ds.rio.write_crs(utm_crs, inplace=True)
        ds_list.append(ds)

        # Print progress
        counter += 1
        print('{:.0f}%'.format((counter / rec_lats.size) * 100), end='\r')

    print('Done')

    # Determine what to output
    if full_output:
        return ds_list, dem_utm
    else:
        return np.array([ds.path_length_difference for ds in ds_list])


def calculate_paths_grid(src_lat, src_lon, radius, spacing, dem_file=None):
    """Calculate paths for a UTM-projected grid surrounding a source location.

    Wrapper around :func:`calculate_paths` for computing path difference grids. See the
    docstring for that function for a description of how DEM data are handled.

    Note:
        Input coordinates are expected to be in the WGS 84 datum. DEM file vertical
        units are expected to be meters.

    Args:
        src_lat (int or float): Source latitude
        src_lon (int or float): Source longitude
        radius (int or float): [m] Desired grid radius, measured from source location
        spacing (int or float): [m] Desired grid spacing
        dem_file (str or None): Path to DEM file (if `None`, then SRTM data are used)

    Returns:
        tuple:  Tuple of the form ``(path_length_differences, dem)`` where
        ``path_length_differences``  is a :class:`~xarray.DataArray` grid of path length
        differences [m], and ``dem`` is a :class:`~xarray.DataArray` containing the
        UTM-projected DEM used to compute the profiles
    """

    # Find UTM CRS of source (see https://gis.stackexchange.com/a/423614)
    utm_crs_list = query_utm_crs_info(
        datum_name='WGS 84',
        area_of_interest=AreaOfInterest(
            west_lon_degree=src_lon,
            south_lat_degree=src_lat,
            east_lon_degree=src_lon,
            north_lat_degree=src_lat,
        ),
    )
    utm_crs = CRS.from_epsg(utm_crs_list[0].code)

    # Get UTM coords for source
    proj = Transformer.from_crs(utm_crs, utm_crs.geodetic_crs)
    src_x, src_y = proj.transform(src_lat, src_lon, direction='INVERSE')

    # Define [gridline-registered] grid of receiver locations [m]
    xlim = (src_x - radius, src_x + radius)
    ylim = (src_y - radius, src_y + radius)

    # Convert gridline registration to pixel registration
    xvec = np.arange(xlim[0] + spacing / 2, xlim[1] + spacing / 2, spacing)
    yvec = np.arange(ylim[0] + spacing / 2, ylim[1] + spacing / 2, spacing)

    # Convert "receiver" UTM coordinates to lat/lon grid
    rec_lat, rec_lon = proj.transform(*np.meshgrid(xvec, yvec))

    # Call calculate_paths()
    tic = time.time()
    ds_list, dem = calculate_paths(
        src_lat=src_lat,
        src_lon=src_lon,
        rec_lat=rec_lat.flatten(),
        rec_lon=rec_lon.flatten(),
        dem_file=dem_file,
        full_output=True,
    )
    toc = time.time()
    print(f'\nElapsed time = {toc - tic:.0f} s')

    # Form a nicely-labeled DataArray from grid of path length differences
    units = dict(units='m')
    path_length_differences = xr.DataArray(
        np.reshape([ds.path_length_difference for ds in ds_list], rec_lat.shape).T,
        coords=[('x', xvec, units), ('y', yvec, units)],
        name='path_length_difference',
        attrs=dict(spacing=spacing, **units),
    ).transpose()
    path_length_differences.rio.write_crs(utm_crs, inplace=True)

    return path_length_differences, dem
