import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

project = 'infresnel'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinxcontrib.apidoc',
]

html_theme = 'sphinx_rtd_theme'
html_show_copyright = False

apidoc_module_dir = '../infresnel'
apidoc_toc_file = False  # Don't create a table of contents file (modules.rst)
apidoc_separate_modules = True  # Give submodules their own page

autodoc_mock_imports = ['numpy', 'pandas', 'pygmt']

napoleon_numpy_docstring = False  # We are using Google docstring style

# These only need to cover the packages we reference from the docstrings
intersphinx_mapping = {
    'numpy': ('https://numpy.org/doc/stable/', None),
    'python': ('https://docs.python.org/3/', None),
    'xarray': ('https://docs.xarray.dev/en/stable/', None),
}
