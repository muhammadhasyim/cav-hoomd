# Configuration file for the Sphinx documentation builder.
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from pathlib import Path

# -- Path setup --------------------------------------------------------------

# Add the project root to the Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / 'src'))
sys.path.insert(0, str(project_root / 'examples'))

# -- Project information -----------------------------------------------------

project = 'Cavity HOOMD'
copyright = '2025, Cavity HOOMD Development Team'
author = 'Cavity HOOMD Development Team'

# The short X.Y version
version = '1.0'
# The full version, including alpha/beta/rc tags
release = '1.0.0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.githubpages',
    'sphinx_autodoc_typehints',
    'myst_parser',
    'sphinx_copybutton',
    'sphinxcontrib.bibtex',
    'nbsphinx',
    'sphinx_design',
    'sphinx_tabs.tabs',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = [
    '_build', 
    'Thumbs.db', 
    '.DS_Store',
    '**.ipynb_checkpoints',
    'examples/.ipynb_checkpoints',
]

# The suffix(es) of source filenames.
source_suffix = {
    '.rst': None,
    '.md': None,
    '.ipynb': None,
}

# The master toctree document.
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.
html_theme = 'furo'

# Theme options are theme-specific and customize the look and feel of a theme
# Following HOOMD-blue color scheme and styling
html_theme_options = {
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
    "top_of_page_button": "edit",
    "light_css_variables": {
        # Primary brand colors (HOOMD-blue inspired)
        "color-brand-primary": "#2980b9",
        "color-brand-content": "#2980b9",
        "color-brand-visited": "#9b59b6",
        
        # Background colors
        "color-background-primary": "#ffffff",
        "color-background-secondary": "#f8f9fa",
        "color-background-hover": "#efeff4",
        "color-background-hover--transparent": "#efeff400",
        "color-background-border": "#e1e4e5",
        
        # Foreground colors
        "color-foreground-primary": "#000000",
        "color-foreground-secondary": "#5a5c63",
        "color-foreground-muted": "#6c6c6c",
        "color-foreground-border": "#878787",
        
        # API documentation colors
        "color-api-background": "#f8f9fa",
        "color-api-background-hover": "#efeff4",
        "color-api-overall": "#2980b9",
        "color-api-name": "#2c3e50",
        "color-api-pre-name": "#34495e",
        "color-api-paren": "#6c6c6c",
        "color-api-keyword": "#8e44ad",
        
        # Admonition colors
        "color-admonition-background": "transparent",
        
        # Inline code colors
        "color-inline-code-background": "#f1f2f3",
        
        # Highlighted colors
        "color-highlighted-background": "#ffd43b",
        "color-highlighted-text": "#000000",
        
        # Guilabel colors
        "color-guilabel-background": "#ddeeff",
        "color-guilabel-border": "#b3d4ff",
        "color-guilabel-text": "#2c3e50",
        
        # Kbd colors
        "color-kbd-background": "#f9f9f9",
        "color-kbd-border": "#d1d5da",
        "color-kbd-text": "#2c3e50",
        
        # Link colors
        "color-link": "#2980b9",
        "color-link-underline": "#2980b9",
        "color-link-hover": "#1abc9c",
        "color-link-hover-underline": "#1abc9c",
    },
    "dark_css_variables": {
        # Primary brand colors (HOOMD-blue dark theme)
        "color-brand-primary": "#3498db",
        "color-brand-content": "#3498db",
        "color-brand-visited": "#bb86fc",
        
        # Background colors
        "color-background-primary": "#131416",
        "color-background-secondary": "#1a1c23",
        "color-background-hover": "#1e2124",
        "color-background-hover--transparent": "#1e212400",
        "color-background-border": "#303335",
        
        # Foreground colors
        "color-foreground-primary": "#ffffffcc",
        "color-foreground-secondary": "#9ca0a5",
        "color-foreground-muted": "#81868d",
        "color-foreground-border": "#666666",
        
        # API documentation colors
        "color-api-background": "#1a1c23",
        "color-api-background-hover": "#1e2124",
        "color-api-overall": "#3498db",
        "color-api-name": "#ecf0f1",
        "color-api-pre-name": "#bdc3c7",
        "color-api-paren": "#95a5a6",
        "color-api-keyword": "#e74c3c",
        
        # Admonition colors
        "color-admonition-background": "transparent",
        
        # Inline code colors
        "color-inline-code-background": "#2c2c2c",
        
        # Highlighted colors
        "color-highlighted-background": "#f39c12",
        "color-highlighted-text": "#000000",
        
        # Guilabel colors
        "color-guilabel-background": "#2c3e50",
        "color-guilabel-border": "#34495e",
        "color-guilabel-text": "#ecf0f1",
        
        # Kbd colors
        "color-kbd-background": "#2c2c2c",
        "color-kbd-border": "#525252",
        "color-kbd-text": "#ecf0f1",
        
        # Link colors
        "color-link": "#3498db",
        "color-link-underline": "#3498db",
        "color-link-hover": "#1abc9c",
        "color-link-hover-underline": "#1abc9c",
    },
    "source_repository": "https://github.com/muhammadhasyim/cav-hoomd/",
    "source_branch": "main",
    "source_directory": "docs/",
    "edit_page_button": True,
    "repository_url": "https://github.com/muhammadhasyim/cav-hoomd/",
    "repository_branch": "main",
    "path_to_docs": "docs",
    "use_edit_page_button": True,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom CSS files
html_css_files = [
    'custom.css',
]

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
# html_logo = '_static/logo.png'

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
# html_favicon = '_static/favicon.ico'

# -- Extension configuration -------------------------------------------------

# -- Options for autodoc extension ------------------------------------------

autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

autodoc_typehints = 'description'
autodoc_typehints_description_target = 'documented'

# -- Options for autosummary extension --------------------------------------

autosummary_generate = True
autosummary_generate_overwrite = True

# -- Options for intersphinx extension --------------------------------------

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'hoomd': ('https://hoomd-blue.readthedocs.io/en/stable/', None),
}

# -- Options for napoleon extension -----------------------------------------

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# -- Options for todo extension ---------------------------------------------

todo_include_todos = True

# -- Options for mathjax extension ------------------------------------------

mathjax3_config = {
    'tex': {
        'inlineMath': [['$', '$'], ['\\(', '\\)']],
        'displayMath': [['$$', '$$'], ['\\[', '\\]']],
        'packages': ['base', 'ams', 'autoload', 'color', 'physics']
    }
}

# -- Options for nbsphinx extension -----------------------------------------

nbsphinx_execute = 'never'  # Don't execute notebooks during build
nbsphinx_allow_errors = True
nbsphinx_timeout = 60

# -- Options for bibtex extension -------------------------------------------

bibtex_bibfiles = ['references.bib']
bibtex_default_style = 'alpha'

# -- Options for copybutton extension ---------------------------------------

copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True
copybutton_line_continuation_character = "\\"

# -- Custom configuration ---------------------------------------------------

# Add custom roles
def setup(app):
    app.add_css_file('custom.css') 