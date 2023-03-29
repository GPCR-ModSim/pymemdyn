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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'pymemdyn'
copyright = '2023, H. Gutierrez de Teran,  X. Bello, M. Esguerra, R. L. van den Broek, R.V. Küpper'
author = 'H. Gutierrez de Teran, X. Bello, M. Esguerra, R. L. van den Broek, R.V. Küpper'

master_doc = 'index'
latex_documents = [
    (master_doc, 'pymemdyn.tex', 'pymemdyn',
        author.replace(', ', ' \\and ').replace(' and ', ' \\and and '),
     'manual'),
]
latex_elements = {'preamble': '\\global\\renewcommand{\\AA}{\\text{\\r{A}}}'}

# The full version, including alpha/beta/rc tags
release = '1.6.1'


# -- General configuration ---------------------------------------------------

import sys

sys.path.append('../')

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
]
autodoc_member_order = 'bysource'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']