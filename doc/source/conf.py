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
# import os
# import sys

from collections import OrderedDict
import MDAnalysis as mda
# import subprocess
import sphinx_rtd_theme
from ipywidgets.embed import DEFAULT_EMBED_REQUIREJS_URL

# -- Project information -----------------------------------------------------

project = 'MDAnalysis User Guide'
copyright = '2019, Lily Wang, Richard J Gowers, Oliver Beckstein'
author = 'Lily Wang, Richard J Gowers, Oliver Beckstein'

# -- Scripts -----------------------------------------------
# Get Travis to regenerate txt tables by re-running scripts
# before deploying docs.
# Turned off for now as 0.21.0 not released yet.
# This allows us to gitignore .txt files as well as
# auto-recreate tables for each deployment.
# Turn off if using sphinx_autobuild as this will autobuild
# to infinity.
#
# subprocess.call('./scripts/generate_all.sh')

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.githubpages',
    'sphinx_sitemap',
    'nbsphinx',
    'sphinx_rtd_theme',
    'IPython.sphinxext.ipython_console_highlighting',
    'IPython.sphinxext.ipython_directive',
    'sphinxcontrib.bibtex',
]

pygments_style = 'default'

todo_include_todos = True

mathjax_path = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML'

site_url = "https://www.mdanalysis.org/UserGuide"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store',
                    '.ipynb_checkpoints', '**/.ipynb_checkpoints', 'scripts', '.*.ipynb']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
# html_theme_path = ['_themes', ]

# styles/fonts to match http://mdanalysis.org (see public/css)
#
# /* MDAnalysis orange: #FF9200 */
# /* MDAnalysis gray: #808080 */
# /* MDAnalysis white: #FFFFFF */
# /* MDAnalysis black: #000000 */

color = {'orange': '#FF9200',
         'gray': '#808080',
         'white': '#FFFFFF',
         'black': '#000000', }

extra_nav_links = OrderedDict()
extra_nav_links['MDAnalysis'] = 'http://mdanalysis.org'
extra_nav_links['docs'] = 'http://docs.mdanalysis.org'
extra_nav_links['wiki'] = 'http://wiki.mdanalysis.org'
extra_nav_links['user discussion group'] = 'http://users.mdanalysis.org'
extra_nav_links['GitHub'] = 'https://github.com/mdanalysis'
extra_nav_links['@mdanalysis'] = 'https://twitter.com/mdanalysis'

html_theme_options = {
    'canonical_url': '',
    'logo_only': True,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'style_nav_header_background': 'white',  # '#e76900', # dark orange
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = "_static/logos/mdanalysis-logo.ico"
html_logo = '_static/logos/user_guide.png'


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_css_files = ['custom.css']  # , 'readable.css']

# Custom sidebar templates, maps document names to template names.
# alabaster sidebars
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',
        'searchbox.html',
    ]
}

# Configuration for intersphinx: refer to the Python standard library
# and other packages used by MDAnalysis
intersphinx_mapping = {'https://docs.python.org/': None,
                       'https://docs.scipy.org/doc/numpy/': None,
                       'https://www.mdanalysis.org/docs/': None,
                       'https://docs.pytest.org/en/latest/': None,
                       }

# nbsphinx
html_js_files = [
    'https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.4/require.min.js',
    DEFAULT_EMBED_REQUIREJS_URL,
]

ipython_warning_is_error = False
nbsphinx_prolog = r"""
.. raw:: html

    <script src='http://cdnjs.cloudflare.com/ajax/libs/require.js/2.1.10/require.min.js'></script>
    <script>require=requirejs;</script>


"""

# substitutions
MDAnalysis_version = '1.0.0'

# rst-epilog implements substitutions
rst_epilog = """
.. |MDAnalysis_version| replace:: {0}
""".format(MDAnalysis_version)
