# import sys
# from pathlib import Path
#
# sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

project = 'infresnel'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinxcontrib.apidoc',
]

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'

html_show_copyright = False

apidoc_module_dir = '../infresnel'

apidoc_toc_file = False  # Don't create a table of contents file (modules.rst)

apidoc_separate_modules = True  # Give submodules their own page

napoleon_numpy_docstring = False  # We are using Google docstring style

intersphinx_mapping = {
    'numpy': ('https://numpy.org/doc/stable/', None),
    'python': ('https://docs.python.org/3/', None),
    'xarray': ('https://docs.xarray.dev/en/stable/', None),
}
