from importlib.metadata import version

from .infresnel import calculate_paths, calculate_paths_grid

__version__ = version('infresnel')
del version
