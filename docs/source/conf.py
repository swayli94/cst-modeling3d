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
import sys
sys.path.insert(0, os.path.abspath('../../'))


# -- Project information -----------------------------------------------------

project = 'CST Modeling'
copyright = '2023, Li Runze'
author = 'Li Runze'

# The full version, including alpha/beta/rc tags
release = 'v0.2.2'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',       # automatically insert docstrings from modules
    'sphinx.ext.mathjax',       # include math, rendered in the browser by MathJax
    'sphinx_rtd_theme',
    'sphinx.ext.napoleon',      # Google' style docstring
    'sphinx.ext.doctest',       # automatically test code snippets in doctest blocks
    'sphinx.ext.intersphinx',   # link between Sphinx documentation of different projects
    'sphinx.ext.todo',          # write "todo" entries that can be shown or hidden on build
    'sphinx.ext.coverage',      # checks for documentation coverage
    'sphinx.ext.viewcode',      # include links to the source code of documented Python objects
    'sphinx.ext.githubpages',   # create .nojekyll file to publish the document on GitHub pages
    'sphinx.ext.napoleon'       # 
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_css_files = ['css/custom.css']

source_encoding = 'utf-8-sig'


# ----------------------------------------------------------------------------
# sphinx_rtd_theme
# https://sphinx-rtd-theme.readthedocs.io/en/stable/configuring.html

html_theme_options = {
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': False,
    'navigation_depth': 5,
    'includehidden': True,
    'titles_only': False
}


# ----------------------------------------------------------------------------

# 激活 图、表、代码块、公式 的自动编号
# 仅针对 有 caption (图例) 标签的对象，该对象的 `numref` 同时生效
numfig = True

# 图例形式
numfig_format = {
    'figure': 'Figure %s',
    'table': 'Table %s',
    'code-block': 'Algorithm %s',
    'section': 'Section %s',
}

# 设置公式编号形式, 如 Eq.10.
math_eqref_format = 'Eq.{number}'  

# 设置所有公式自动编号
# 否则需自己标注 :label:
math_number_all = False

# 设置公式编号包含的章节层级
math_numfig = True
numfig_secnum_depth = 2

