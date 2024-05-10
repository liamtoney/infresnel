import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pygmt
import xarray as xr
from joblib import Parallel, delayed
from pyproj import Transformer
from rasterio.enums import Resampling
from scipy.interpolate import RectBivariateSpline
from tqdm.auto import tqdm

from ._georeference import (
    _check_valid_elevation_for_coords,
    _estimate_utm_crs,
    _export_geotiff,
)
from ._path import (
    _direct_path,
    _horizontal_distance,
    _norm,
    _path_length,
    _shortest_diffracted_path,
)


def calculate_paths(
    src_lat,
    src_lon,
    rec_lat,
    rec_lon,
    dem_file=None,
    full_output=False,
    return_dem=False,
    n_jobs=1,
):
    """Calculate elevation profiles, direct paths, and shortest diffracted paths.

    Paths are calculated for a given DEM (either a user-supplied `dem_file` or
    automatically downloaded 1 arc-second SRTM data) and an arbitrary number of
    source-receiver pairs. By default, the function returns only the direct and shortest
    diffracted path lengths. If `full_output` is set to `True`, then complete path
    information (lengths, coordinates, etc.) is returned.

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
        full_output (bool): Toggle outputting full profile/path information vs. just
            direct and shortest diffracted path lengths
        return_dem (bool): Toggle additionally returning the UTM-projected DEM used to
            compute the profiles
        n_jobs (int): Number of parallel jobs to run (default is 1, which means no
            parallelization) — this argument is passed on to :class:`joblib.Parallel`

    Returns:
        If `full_output` is `False` — an :class:`~numpy.ndarray` with shape ``(2,
        n_receivers)`` containing the direct path lengths [m] (first row) and shortest
        diffracted path lengths [m] (second row) for each source-receiver pair.

        If `full_output` is `True` — a list of :class:`~xarray.Dataset` objects, one per
        source-receiver pair, containing full profile and path information

        For either of the above cases, if `return_dem` is `True`, then the function
        returns a tuple of the form ``(output_array, dem)`` where ``output_array``  is
        the output described above, and ``dem`` is a :class:`~xarray.DataArray`
        containing the UTM-projected DEM used to compute the profiles
    """

    # Type checks
    message = 'src_lat and src_lon must both be scalars!'
    assert np.isscalar(src_lat) and np.isscalar(src_lon), message

    # Type conversion, so we can iterate
    rec_lats = np.atleast_1d(rec_lat)
    rec_lons = np.atleast_1d(rec_lon)

    # Define number of paths
    n_paths = rec_lats.size

    print('Loading and projecting DEM...')
    if dem_file is not None:
        # Load user-provided DEM, first checking if it exists
        dem_file = Path(str(dem_file)).expanduser().resolve()
        assert dem_file.is_file(), 'dem_file does not exist!'
        dem = xr.open_dataarray(dem_file)
        # User-provided DEM may not fully encompass source and receivers!
        sufficient_extent_guaranteed = False
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
        # PyGMT fills the "nodata" area with zeros, which is same as water. So we
        # convert these areas to NaN here. See
        # https://en.wikipedia.org/wiki/Shuttle_Radar_Topography_Mission for discussion
        # of the data bounds of the SRTM data.
        dem = dem.where((dem.lat < 60) & (dem.lat > -56))
        dem.rio.write_nodata(np.nan, inplace=True)
        # We know that this DEM fully encompasses source and receivers (by design)!
        sufficient_extent_guaranteed = True

    # Clean DEM before going further
    dem = dem.squeeze(drop=True).rename('elevation')

    # Project DEM to UTM, further relabeling
    utm_crs = _estimate_utm_crs(src_lat, src_lon, datum_name='WGS 84')
    dem_utm = dem.rio.reproject(utm_crs, resampling=Resampling.cubic_spline)
    units = dict(units='m')
    dem_utm.attrs = units
    for coordinate in 'x', 'y':
        dem_utm[coordinate].attrs = units
    print('Done\n')

    # Evaluate presence of NaN values in DEM; determine if we need to run "costly check"
    if dem.isnull().all():
        raise ValueError('DEM is entirely NaN values! Exiting.')
    elif dem.isnull().any():
        warnings.warn(
            f'{dem.isnull().values.sum() / dem.size:.1%} of DEM is NaN!', stacklevel=2
        )
        sys.stderr.flush()
        check_for_valid_elevations = True  # Since we have NaNs in DEM, we should check
    else:  # DEM values are all valid
        # Even if the DEM is fully valid, for user-supplied DEMs the extent might not
        # cover all sources and receivers, so in that case we still should check
        check_for_valid_elevations = not sufficient_extent_guaranteed

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

    # Check that source and receiver(s) all have valid elevation values in DEM (this is
    # the "costly check" mentioned above). Mainly relevant for user-supplied DEMs... but
    # also must be run if the PyGMT-supplied DEM is not fully within SRTM range (kind of
    # unlikely edge case). For most PyGMT-supplied DEMs, this check will not end up
    # being run — which is good, since it can be SLOW.
    if check_for_valid_elevations:
        print('Checking that DEM contains source and receivers...')
        if not _check_valid_elevation_for_coords(
            dem_utm, mean_resolution, src_x, src_y
        ):
            raise ValueError('Source is not in DEM! Exiting.')

        # KEY: Parallel computation of valid receiver paths
        compute_paths = np.array(
            Parallel(n_jobs=n_jobs)(
                delayed(_check_valid_elevation_for_coords)(
                    dem_utm, mean_resolution, rec_x, rec_y
                )
                for rec_x, rec_y in zip(rec_xs, rec_ys)
            )
        )

        n_invalid_paths = (~compute_paths).sum()
        if n_invalid_paths > 0:
            print(
                f'Done — {n_invalid_paths} invalid path{"" if n_invalid_paths == 1 else "s"} will be set to NaN\n'
            )
        else:
            print('Done\n')
    else:
        compute_paths = np.full(n_paths, True)  # Compute all paths

    # Fit bivariate spline to DEM (slow for very high resolution DEMs!)
    print('Fitting spline to DEM...')
    x = dem_utm.x
    y = dem_utm.y
    z = dem_utm.fillna(dem_utm.median())  # Can't have NaNs in z
    if not pd.Series(x).is_monotonic_increasing:
        x = x[::-1]
        z = np.fliplr(z)
    if not pd.Series(y).is_monotonic_increasing:
        y = y[::-1]
        z = np.flipud(z)
    spline = RectBivariateSpline(x=x, y=y, z=z.T)  # x and y are monotonic increasing
    print('Done\n')

    # Define fuction for calculating a single path
    def _calculate_single_path(rec_x, rec_y, compute_path):
        # If the DEM points were valid, compute the path
        if compute_path:
            # Determine # of points in profile
            dist = _norm(src_x - rec_x, src_y - rec_y)
            n = max(int(np.ceil(dist / target_spacing)), 2)  # Ensure at least 2 points!

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
            diff_path = _shortest_diffracted_path(
                profile.distance.values, profile.values
            )
            diff_path_len = _path_length(profile.distance.values, diff_path)

        # Just populate everything with NaNs
        else:
            profile = direct_path = diff_path = [np.nan]
            direct_path_len = diff_path_len = np.nan

        # Choose output format
        if full_output:
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
                attrs=dict(
                    path_length_difference=diff_path_len - direct_path_len, **units
                ),
            )
            ds.rio.write_crs(utm_crs, inplace=True)
            output = ds  # KEY: Return a Dataset
        else:
            # Just include the path lengths
            output = direct_path_len, diff_path_len  # KEY: Return a tuple

        return output

    print(f'Computing {n_paths} path{"" if n_paths == 1 else "s"}...')

    # Only create the progress bar if we have more than 1 path
    iterable = zip(rec_xs, rec_ys, compute_paths)
    if n_paths > 1:
        iterable = tqdm(
            iterable,
            total=n_paths,
            bar_format='{percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} paths ',
        )

    # KEY: Parallel path calculations over all receivers (= source-receiver pairs)
    output_array = Parallel(n_jobs=n_jobs)(
        delayed(_calculate_single_path)(rec_x, rec_y, compute_path)
        for rec_x, rec_y, compute_path in iterable
    )

    print('Done')

    # Determine what to output
    if not full_output:
        output_array = np.array(output_array).T  # Convert list of tuples to array
    if return_dem:
        return output_array, dem_utm
    else:
        return output_array


def calculate_paths_grid(
    src_lat,
    src_lon,
    x_radius,
    y_radius,
    spacing,
    dem_file=None,
    output_file=None,
    n_jobs=1,
):
    """Calculate paths for a UTM-projected grid surrounding a source location.

    Wrapper around :func:`calculate_paths` for computing path difference grids. See the
    docstring for that function for a description of how DEM data are handled.

    Note:
        Input coordinates are expected to be in the WGS 84 datum. DEM file vertical
        units are expected to be meters.

    Args:
        src_lat (int or float): Source latitude
        src_lon (int or float): Source longitude
        x_radius (int, float, list, or tuple): [m] Desired grid radius in
            :math:`x`-direction, measured from source location (specify a two-element
            array for different west and east extents)
        y_radius (int, float, list, or tuple): [m] Desired grid radius in
            :math:`y`-direction, measured from source location (specify a two-element
            array for different south and north extents)
        spacing (int or float): [m] Desired grid spacing
        dem_file (str or None): Path to DEM file (if `None`, then SRTM data are used)
        output_file (str or None): If a string filepath is provided, then an RGB GeoTIFF
            file containing the colormapped grid of path length difference values is
            exported to this filepath (no export if `None`)
        n_jobs (int): Number of parallel jobs to run (default is 1, which means no
            parallelization) — this argument is passed on to :class:`joblib.Parallel`

    Returns:
        tuple:  Tuple of the form ``(path_length_differences, dem)`` where
        ``path_length_differences``  is a :class:`~xarray.DataArray` grid of path length
        differences [m], and ``dem`` is a :class:`~xarray.DataArray` containing the
        UTM-projected DEM used to compute the profiles
    """

    # Find UTM CRS of source
    utm_crs = _estimate_utm_crs(src_lat, src_lon, datum_name='WGS 84')

    # Get UTM coords for source
    proj = Transformer.from_crs(utm_crs, utm_crs.geodetic_crs)
    src_x, src_y = proj.transform(src_lat, src_lon, direction='INVERSE')

    # Function for pre-processing radii arguments
    def _process_radius(radius):
        radius = np.atleast_1d(radius)
        if radius.size == 1:
            radius = radius.repeat(2)
        elif radius.size == 2:
            pass  # We already have a two-element array!
        else:
            raise ValueError('x_radius and y_radius each take only one or two values!')
        return radius

    # Define [gridline-registered] grid of receiver locations [m]
    x_radius = _process_radius(x_radius)
    y_radius = _process_radius(y_radius)
    xlim = (src_x - x_radius[0], src_x + x_radius[1])
    ylim = (src_y - y_radius[0], src_y + y_radius[1])

    # Convert gridline registration to pixel registration
    xvec = np.arange(xlim[0] + spacing / 2, xlim[1] + spacing / 2, spacing)
    yvec = np.arange(ylim[0] + spacing / 2, ylim[1] + spacing / 2, spacing)
    yvec = yvec[::-1]  # This places the origin at top-left

    # Convert "receiver" UTM coordinates to lat/lon grid
    rec_lat, rec_lon = proj.transform(*np.meshgrid(xvec, yvec))

    # Call calculate_paths()
    tic = time.time()
    (direct_path_lens, diff_path_lens), dem = calculate_paths(
        src_lat=src_lat,
        src_lon=src_lon,
        rec_lat=rec_lat.flatten(),
        rec_lon=rec_lon.flatten(),
        dem_file=dem_file,
        return_dem=True,
        n_jobs=n_jobs,
    )
    toc = time.time()
    print(f'\nElapsed time = {toc - tic:.0f} s')

    # Form a nicely-labeled DataArray grid of path length differences
    units = dict(units='m')
    path_length_differences = xr.DataArray(
        np.reshape(diff_path_lens - direct_path_lens, rec_lat.shape).T,
        coords=[('x', xvec, units), ('y', yvec, units)],
        name='path_length_difference',
        attrs=dict(spacing=spacing, **units),
    ).transpose()
    path_length_differences.rio.write_crs(utm_crs, inplace=True)

    # Export GeoTIFF if requested
    if output_file is not None:
        print()
        _export_geotiff(path_length_differences, filename=output_file)

    return path_length_differences, dem
