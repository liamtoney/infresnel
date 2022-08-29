"""
Reproduce Fig. 2B in Fee et al. (2021) (https://doi.org/10.3389/feart.2021.620813) using
path difference time delays.

TODO: Make this into a proper testing file for use with pytest...
"""

import time

# We need Matplotlib, which is an optional dependency - so we install here if needed
# fmt: off
try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    import subprocess
    subprocess.run(['pip', 'install', 'matplotlib'])
    import matplotlib.pyplot as plt
# fmt: on
import numpy as np
import pandas as pd
from pyproj import CRS, Transformer

from infresnel import calculate_paths

CELERITY = 343.5  # [m/s] From Fee et al. (2021)

# Full path to DEM_Union_UAV_161116_sm101.tif (121 MB GeoTIFF, get from Liam or David)
DEM_FILE = '/Users/ldtoney/work/yasur_ml/data/DEM_Union_UAV_161116_sm101.tif'

# Load YIF6 coordinates (we are using YIF6 as a source!)
SRC_URL = 'http://service.iris.edu/fdsnws/station/1/query?net=3E&sta=YIF6&format=geocsv'
src_ser = pd.read_csv(SRC_URL, sep='|', comment='#').squeeze()
src_lat, src_lon = src_ser.Latitude, src_ser.Longitude

# Hard-coded for Yasur, UTM zone 59S
YASUR_CRS = CRS(32759)

# Define [gridline-registered] grid of receiver locations
SPACING = 100  # [m] 67 s
# SPACING = 50  # [m] 249 s
# SPACING = 25  # [m] 968 s
# SPACING = 10  # [m] 6201 s
XLIM = (336000, 338200)  # [m] in YASUR_CRS
YLIM = (7839200, 7840600)  # [m] in YASUR_CRS

# Converting gridline registration to pixel registration below
xvec = np.arange(XLIM[0] + SPACING / 2, XLIM[1] + SPACING / 2, SPACING)
yvec = np.arange(YLIM[0] + SPACING / 2, YLIM[1] + SPACING / 2, SPACING)
xgrid, ygrid = np.meshgrid(xvec, yvec)

# Convert "receiver" UTM coordinates to lat/lon
proj = Transformer.from_crs(YASUR_CRS, YASUR_CRS.geodetic_crs)
rec_lat, rec_lon = proj.transform(xgrid, ygrid)

#%% Calculate travel time delay grid

tic = time.time()

path_length_differences = calculate_paths(
    src_lat=src_lat,
    src_lon=src_lon,
    rec_lat=rec_lat.flatten(),
    rec_lon=rec_lon.flatten(),
    dem_file=DEM_FILE,
)

toc = time.time()
print(f'\nElapsed time = {toc - tic:.0f} s')

travel_time_delays = np.reshape(path_length_differences / CELERITY, rec_lat.shape)

#%% Plot travel time delay grid

# Convert from pixel registration back to gridline registration for plotting
xvec_plot = np.hstack([xvec, xvec[-1] + SPACING]) - SPACING / 2
yvec_plot = np.hstack([yvec, yvec[-1] + SPACING]) - SPACING / 2

# Make grid corner (0, 0) by applying this simple transform function
transform = lambda x, y: (x - xvec_plot.min(), y - yvec_plot.min())

fig, ax = plt.subplots(figsize=(8, 4))
im = ax.pcolormesh(
    *transform(xvec_plot, yvec_plot),
    travel_time_delays,
    cmap='magma_r',
    shading='flat',
    # vmax=0.42,  # Match Fig. 2B in Fee et al. (2021)
)
src_utm = proj.transform(src_lat, src_lon, direction='INVERSE')
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
ax.set_title(f'{SPACING} m spacing — {toc - tic:.0f} s')
ax.grid(linestyle=':', alpha=0.5)
ax.set_aspect('equal')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
fig.colorbar(im, label='Time (s)')
fig.show()
