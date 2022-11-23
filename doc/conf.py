# import sys
# from pathlib import Path
#
# sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

project = 'infresnel'
copyright = '2022, Liam Toney'
author = 'Liam Toney'

extensions = ['sphinx.ext.autodoc', 'sphinxcontrib.apidoc', 'sphinx.ext.napoleon']

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'

apidoc_module_dir = '../infresnel'

apidoc_toc_file = False  # Don't create a table of contents file (modules.rst)

napoleon_numpy_docstring = False  # We are using Google docstring style
