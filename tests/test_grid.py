"""
Reproduce Fig. 2B in Fee et al. (2021) (https://doi.org/10.3389/feart.2021.620813) using
path difference time delays.

TODO: Make this into a proper testing file for use with pytest...
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyproj import CRS, Transformer

from infresnel import calculate_paths_grid

CELERITY = 343.5  # [m/s] From Fee et al. (2021)

# Full path to DEM_Union_UAV_161116_sm101.tif (121 MB GeoTIFF, get from Liam or David)
DEM_FILE = '/Users/ldtoney/work/yasur_ml/data/DEM_Union_UAV_161116_sm101.tif'

# Load YIF6 coordinates (we are using YIF6 as a source!)
SRC_URL = 'http://service.iris.edu/fdsnws/station/1/query?net=3E&sta=YIF6&format=geocsv'
src_ser = pd.read_csv(SRC_URL, sep='|', comment='#').squeeze()
src_lat, src_lon = src_ser.Latitude, src_ser.Longitude

# Define [gridline-registered] grid of receiver locations
SPACING = 100  # [m]
X_RADIUS = (450, 1750)  # [m] Hardcoded to match Fig. 2B
Y_RADIUS = (1200, 200)  # [m] "                        "

#%% Calculate travel time delay grid

path_length_differences, dem = calculate_paths_grid(
    src_lat=src_lat,
    src_lon=src_lon,
    x_radius=X_RADIUS,
    y_radius=Y_RADIUS,
    spacing=SPACING,
    dem_file=DEM_FILE,
)

travel_time_delays = path_length_differences.data / CELERITY

#%% Plot travel time delay grid

# Convert from pixel registration back to gridline registration for plotting
xvec, yvec = path_length_differences.x, path_length_differences.y
spacing = path_length_differences.spacing
xvec_plot = np.hstack([xvec, xvec[-1] + spacing]) - spacing / 2
yvec_plot = np.hstack([yvec, yvec[-1] - spacing]) + spacing / 2  # Note sign change!

# Make grid corner (0, 0) by applying this simple transform function
transform = lambda x, y: (x - xvec_plot.min(), y - yvec_plot.min())

fig, ax = plt.subplots(figsize=(8, 4))
hs = dem.copy()
hs.data = matplotlib.colors.LightSource().hillshade(
    dem.data,
    dx=abs(dem.x.diff('x').mean().values),
    dy=abs(dem.y.diff('y').mean().values),
)
hs['x'], hs['y'] = transform(hs.x, hs.y)
hs.plot.imshow(ax=ax, cmap='Greys_r', add_colorbar=False)
im = ax.pcolormesh(
    *transform(xvec_plot, yvec_plot),
    travel_time_delays,
    cmap='magma_r',
    shading='flat',
    # vmax=0.42,  # Match Fig. 2B in Fee et al. (2021)
    alpha=0.8,
)
utm_crs = CRS(path_length_differences.rio.crs)
proj = Transformer.from_crs(utm_crs.geodetic_crs, utm_crs)
src_utm = proj.transform(src_lat, src_lon)
ax.scatter(*transform(*src_utm), s=180, marker='v', color='black', lw=0, zorder=2)
ax.text(
    *transform(src_utm[0], src_utm[1] + 5),  # [m] Shifting text here
    src_ser.Station[-1],  # Just use station # here
    color='white',
    ha='center',
    va='center',
    fontsize=8,
    weight='bold',
    zorder=2,
)
ax.set_title(f'{SPACING} m spacing')
ax.set_xlim(0, xvec_plot.max() - xvec_plot.min())
ax.set_ylim(0, yvec_plot.max() - yvec_plot.min())
ax.grid(linestyle=':', alpha=0.5)
ax.set_aspect('equal')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
fig.colorbar(im, label='Time (s)')
fig.show()
