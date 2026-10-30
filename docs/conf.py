# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# The package is imported from the repository rather than installed, so the
# docs build doesn't need the GUI dependencies.
import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'gateway_mt'
copyright = '2019, Eric Roots'
author = 'Eric Roots'

# The full version, including alpha/beta/rc tags, read from gateway_mt/__init__.py
import re
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gateway_mt', '__init__.py')) as _f:
    release = re.search(r"__version__ = '([^']+)'", _f.read()).group(1)


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx.ext.intersphinx',
              'sphinx.ext.viewcode']

# API pages: list members in source order, and read NumPy-style docstrings
autodoc_member_order = 'bysource'
autodoc_default_options = {'members': True, 'show-inheritance': True}
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_preprocess_types = True     # {'a', 'b'} value sets render as literals
napoleon_use_ivar = True             # class 'Attributes' sections don't clash with properties

intersphinx_mapping = {'python': ('https://docs.python.org/3', None),
                       'numpy': ('https://numpy.org/doc/stable', None),
                       'scipy': ('https://docs.scipy.org/doc/scipy', None),
                       'matplotlib': ('https://matplotlib.org/stable', None)}

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build']
source_suffix = '.rst'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
