# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import subprocess
import sys
# sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../../'))


# -- Project information -----------------------------------------------------

project = 'Materials Learning Algorithms Data Acquisition (MALADA)'
copyright = 'noneyet'

author = 'L. Fiedler'

# The version info for the project
tag = subprocess.run(['git', 'describe', '--tags'], capture_output=True,
                        text=True)
version = tag.stdout.strip()


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'myst_parser',
    'sphinx_markdown_tables',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
#    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
]

napoleon_google_docstring = False
napoleon_numpy_docstring = True

autodoc_mock_imports = [
    'ase',
    'optuna',
    'mpmath',
    'torch',
    'numpy',
    'scipy',
    'oapackage',
    'matplotlib',
    'horovod',
    'lammps',
    'total_energy',
    'asap3'
]

autodoc_member_order = 'groupwise'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

html_context = {
  'display_github': True,
  'github_repo': 'mala-project/malada',
  'github_version': 'develop',
  'conf_py_path': '/docs/source/',
}

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#html_logo = "./img/logos/mala_vertical.png"

# The name of an image file (relative to this directory) to use as a favicon of
# the docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = "./img/logos/mala_favicon.png"

# The suffix of source file names.
source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}
