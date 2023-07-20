from pathlib import Path

import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from pyproj import CRS
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info

# GeoTIFF datatype must be "byte" for Google Earth
BYTE = np.ubyte


# Find UTM CRS of a (lat, lon) point (see https://gis.stackexchange.com/a/423614)
def _estimate_utm_crs(lat, lon, datum_name='WGS 84'):
    utm_crs_list = query_utm_crs_info(
        datum_name=datum_name,
        area_of_interest=AreaOfInterest(
            west_lon_degree=lon,
            south_lat_degree=lat,
            east_lon_degree=lon,
            north_lat_degree=lat,
        ),
    )
    return CRS.from_epsg(utm_crs_list[0].code)  # Taking first entry of list here!


# Export a DataArray grid as an RGB GeoTIFF
def _export_geotiff(grid, filename, cmap=cc.m_fire_r):
    # Process filepath
    if not (filename.endswith('.tif') or filename.endswith('.tiff')):
        filename += '.tiff'
    filename = Path(filename).expanduser().resolve()

    # Map grid values to full colormap (no clipping!), drop alpha, and convert datatype
    color_data = cmap(plt.Normalize()(grid.data))[:, :, :-1]
    color_data = (color_data * np.iinfo(BYTE).max).astype(BYTE)

    # Form RGB DataArray
    da = xr.DataArray(color_data, coords=dict(y=grid.y, x=grid.x, band=range(3)))
    da = da.transpose('band', 'y', 'x')  # This order is required for export

    # Add CRS information and write to file
    da.rio.write_crs(grid.rio.crs, inplace=True)
    da.rio.to_raster(filename, dtype=BYTE, tags=dict(AREA_OR_POINT='Point'))
    print(f'GeoTIFF saved to {filename}')


# Check that a pair of UTM coordinates has a valid (not NaN!) elevation value in a UTM DEM
def _check_valid_elevation_for_coords(dem, tolerance, x, y):
    try:
        # We're *inside* the DEM extent (non-rectangular DEM shapes still might not have
        # a valid elevation value, though!)
        queried_elevation = dem.sel(
            x=x, y=y, method='nearest', tolerance=tolerance
        ).data
    except KeyError:
        # We're *outside* the DEM extent
        queried_elevation = np.nan
    return not np.isnan(queried_elevation)  # Boolean is `True` for valid elevations!
