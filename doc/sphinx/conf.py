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
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'MDAnalysis User Guide'
copyright = '2019, Lily Wang, Richard J Gowers, Oliver Beckstein'
author = 'Lily Wang, Richard J Gowers, Oliver Beckstein'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx_sitemap',
]

mathjax_path = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML'

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
html_theme = 'alabaster'

# styles/fonts to match http://mdanalysis.org (see public/css)
#
# /* MDAnalysis orange: #FF9200 */
# /* MDAnalysis gray: #808080 */
# /* MDAnalysis white: #FFFFFF */
# /* MDAnalysis black: #000000 */

color = {'orange': '#FF9200',
         'gray': '#808080',
         'white': '#FFFFFF',
         'black': '#000000',}

from collections import OrderedDict
extra_nav_links = OrderedDict()
extra_nav_links['MDAnalysis'] = 'http://mdanalysis.org'
extra_nav_links['docs'] = 'http://docs.mdanalysis.org'
extra_nav_links['wiki'] = 'http://wiki.mdanalysis.org'
extra_nav_links['user discussion group'] = 'http://users.mdanalysis.org'
extra_nav_links['GitHub'] = 'https://github.com/mdanalysis'
extra_nav_links['@mdanalysis'] = 'https://twitter.com/mdanalysis'

html_theme_options = {
    'logo' : "logos/mdanalysistutorial-logo-200x150.png",
    'github_user': 'MDAnalysis',
    'github_repo': 'mdanalysis',
    'github_button': False,
    # 'github_type': 'star',
    'github_banner': True,
    'show_related': True,
    'fixed_sidebar': False,
    'sidebar_includehidden': True,
    'sidebar_collapse': True,
    # style
    'link': color['orange'],
    'link_hover': color['orange'],
    'gray_1': color['gray'],
    'narrow_sidebar_bg': color['gray'],
    'narrow_sidebar_fg': color['white'],
    # typography
    'font_family': "'PT Sans', Helvetica, Arial, 'sans-serif'",
    'head_font_family': "",
    'code_font_family': "Menlo, Monaco, 'Courier New', monospace",
    'caption_font_size': "smaller",
    # external links
    'extra_nav_links': extra_nav_links,
}

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = "_static/logos/mdanalysis-logo.ico"


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

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
