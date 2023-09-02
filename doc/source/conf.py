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

import datetime
import MDAnalysis as mda
import subprocess
from ipywidgets.embed import DEFAULT_EMBED_REQUIREJS_URL

# -- Project information -----------------------------------------------------

project = "MDAnalysis User Guide"


def sort_authors(filename):
    """Generate sorted list of authors from AUTHORS"""
    authors = []
    with open(filename, "r") as f:
        contents = f.read()
    lines = contents.split("Chronological list of authors")[1].split("\n")[2:]
    lines = [x for x in lines if x]
    for line in lines:
        line = line.strip()
        if line[:2] == "- ":
            authors.append(line[2:].strip())

    # remove original authors
    original = ["Lily Wang", "Richard J. Gowers", "Oliver Beckstein"]
    for name in original:
        authors.remove(name)

    # sort on last name
    authors.sort(key=lambda name: name.split()[-1])
    authors = original[:1] + authors + original[-2:]
    return authors


author_list = sort_authors("AUTHORS")
author = ", ".join(author_list[:-1]) + ", and " + author_list[-1]
now = datetime.datetime.now()
copyright = "2019-{}, {}.".format(now.year, author)

# -- Scripts -----------------------------------------------
# Get Travis to regenerate txt tables by re-running scripts
# before deploying docs.
# This allows us to gitignore .txt files as well as
# auto-recreate tables for each deployment.
# Turn off if using sphinx_autobuild as this will autobuild
# to infinity.
# TURNED OFF FOR NOW until we sort out how versioned docs
# work and how we will install MDAnalysis + dependencies
#
subprocess.call("./scripts/generate_all.sh")

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.githubpages",
    "sphinx_sitemap",
    "nbsphinx",
    "sphinx_rtd_theme",
    "IPython.sphinxext.ipython_console_highlighting",
    "IPython.sphinxext.ipython_directive",
    "sphinxcontrib.bibtex",
    "matplotlib.sphinxext.plot_directive",
]

bibtex_bibfiles = ["references.bib"]

pygments_style = "default"

todo_include_todos = True

mathjax_path = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

html_baseurl = "https://userguide.mdanalysis.org/"
sitemap_url_scheme = "{link}"
html_use_opensearch = "https://userguide.mdanalysis.org"


# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    ".ipynb_checkpoints",
    "**/.ipynb_checkpoints",
    "scripts",
    "**/.*.ipynb",
    ".*",
]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'mdanalysis_sphinx_theme'


# styles/fonts to match http://mdanalysis.org (see public/css)
#
# /* MDAnalysis orange: #FF9200 */
# /* MDAnalysis gray: #808080 */
# /* MDAnalysis white: #FFFFFF */
# /* MDAnalysis black: #000000 */

color = {
    "orange": "#FF9200",
    "gray": "#808080",
    "white": "#FFFFFF",
    "black": "#000000",
}

extra_nav_links = {}
extra_nav_links['MDAnalysis'] = 'http://mdanalysis.org'
extra_nav_links['docs'] = 'http://docs.mdanalysis.org'
extra_nav_links['wiki'] = 'http://wiki.mdanalysis.org'
extra_nav_links['user discussion group'] = 'http://users.mdanalysis.org'
extra_nav_links['GitHub'] = 'https://github.com/mdanalysis'
extra_nav_links['@mdanalysis'] = 'https://twitter.com/mdanalysis'

html_theme_options = {
    "canonical_url": "",
    "logo_only": True,
    "display_version": True,
    "prev_next_buttons_location": "bottom",
    "style_external_links": False,
    "style_nav_header_background": "white",  # '#e76900', # dark orange
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
    'mda_official': True,
}

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = "_static/logos/mdanalysis-logo.ico"
html_logo = "_static/logos/user_guide.png"

html_context = {"versions_json_url": "https://userguide.mdanalysis.org/versions.json"}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
# html_css_files = ['custom.css']  # , 'readable.css']

# Custom sidebar templates, maps document names to template names.
# alabaster sidebars
html_sidebars = {
    "**": [
        "about.html",
        "navigation.html",
        "relations.html",
        "searchbox.html",
    ]
}

# Configuration for intersphinx: refer to the Python standard library
# and other packages used by MDAnalysis
mda_version = mda.__version__
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "mdanalysis": (f"https://docs.mdanalysis.org/{mda_version}/", None),
    "pytest": ("https://docs.pytest.org/en/latest/", None),
    "chemfiles": ("https://chemfiles.org/chemfiles.py/latest/", None),
    "parmed": ("https://parmed.github.io/ParmEd/html/", None),
}

# nbsphinx
html_js_files = [
    # 'https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.4/require.min.js',
    # DEFAULT_EMBED_REQUIREJS_URL,
]

ipython_warning_is_error = False
nbsphinx_prolog = r"""
.. raw:: html

    <script src='http://cdnjs.cloudflare.com/ajax/libs/require.js/2.1.10/require.min.js'></script>
    <script>require=requirejs;</script>

"""

# substitutions
MDAnalysis_version = version = mda.__version__

# rst-epilog implements substitutions
rst_epilog = f"""
.. |MDAnalysis_version| replace:: {MDAnalysis_version}
"""
