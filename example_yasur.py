from pathlib import Path
from urllib.request import urlretrieve

import pandas as pd

from infresnel import calculate_paths

# Load Yasur subcrater locations
SRC_URL = 'https://raw.githubusercontent.com/liamtoney/yasur_ml/main/yasur_subcrater_locs.json'
src_df = pd.read_json(SRC_URL)

# Load Yasur station coordinates
REC_URL = 'http://service.iris.edu/fdsnws/station/1/query?net=3E&sta=YIF?&format=geocsv'
rec_df = pd.read_csv(REC_URL, sep='|', comment='#')
rec_df = rec_df[rec_df.Station != 'YIF6']  # Don't include YIF6

# Download Yasur DEM, if it doesn't already exist (~200 MB GeoTIFF)
DEM_URL = 'https://opentopography.s3.sdsc.edu/dataspace/OTDS.072019.4326.1/raster/DEM_WGS84.tif'
dem_file = DEM_URL.split('/')[-1]
if not Path(dem_file).exists():
    print('Downloading DEM...')
    urlretrieve(DEM_URL, dem_file)
    print('Done\n')

# Call function
ds_list, dem = calculate_paths(
    src_lat=src_df.S[1],
    src_lon=src_df.S[0],
    rec_lat=rec_df.Latitude,
    rec_lon=rec_df.Longitude,
    dem_file=dem_file,
    full_output=True,
)

#%% Create example figures (optional)

SAVE_EXAMPLE_FIGURES = False

# We need Matplotlib, which is an optional dependency - so we install here if needed
# fmt: off
try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    import subprocess
    subprocess.run(['pip', 'install', 'matplotlib'])
    import matplotlib.pyplot as plt
# fmt: on

# Reset everything to defaults; use smaller font size
plt.rcParams.update(plt.rcParamsDefault)
plt.rc('font', size=9)

# Plot DEM with source-receiver paths
fig, ax = plt.subplots()
dem.plot.imshow(
    ax=ax, cmap='Greys_r', center=False, cbar_kwargs=dict(label='Elevation (m)')
)
for ds, station in zip(ds_list, rec_df.Station):
    ax.plot(ds.x, ds.y, solid_capstyle='round', label=station)
ax.scatter(ds.x[0], ds.y[0], c='white', ec='black', zorder=2, label='Source')
ax.ticklabel_format(style='plain')
ax.set_aspect('equal')
ax.set_xlabel('UTM easting (m)')
ax.set_ylabel('UTM northing (m)')
ax.legend(loc='lower right', frameon=False)
fig.tight_layout()
fig.show()
if SAVE_EXAMPLE_FIGURES:
    fig.savefig('example_figures/yasur_dem_paths.png', bbox_inches='tight', dpi=300)

# Plot comparison of elevation profiles, direct paths, and shortest diffracted paths
fig, axes = plt.subplots(nrows=3, sharex=True, sharey=True)
for ax, var_name in zip(axes, ds.data_vars):
    for ds, station in zip(ds_list, rec_df.Station):
        ax.plot(ds.distance, ds[var_name], solid_capstyle='round', label=station)
    ax.scatter(
        ds.distance[0], ds.elevation[0], c='white', ec='black', zorder=2, label='Source'
    )
    ax.set_title(var_name, fontsize='medium', fontname='monospace')
    ax.set_aspect('equal')
axes[-1].set_xlabel('Horizontal distance (m)')
axes[1].set_ylabel('Elevation (m)')
axes[1].legend(loc='center left', frameon=False, bbox_to_anchor=(1.05, 0.5))
fig.tight_layout()
fig.show()
if SAVE_EXAMPLE_FIGURES:
    fig.savefig(
        'example_figures/yasur_path_comparison.png', bbox_inches='tight', dpi=300
    )
