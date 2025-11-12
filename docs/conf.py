# Configuration file for the Sphinx documentation builder.
#
# SIMPLIFIED ARCHITECTURE:
# This configuration uses a streamlined approach after project restructuring:
# 1. HOOMD is installed via conda (see environment.yml) 
# 2. Only C++ extensions (_cavitymd, _bussi_reservoir) are mocked
# 3. Plugins are imported directly from src/ (no intermediate hoomd/ folder)
# 4. Real docstrings and signatures are preserved
# 
# This eliminates the complex mocking system previously required.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from pathlib import Path

# Add the project paths to Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / 'src'))
sys.path.insert(0, str(project_root / 'examples'))

# Check if we need to mock imports (either on RTD or any CI environment without HOOMD)
on_rtd = (os.environ.get('READTHEDOCS') == 'True' or 
          os.environ.get('RTD_ENV_NAME') is not None or
          'readthedocs' in os.environ.get('HOSTNAME', '').lower())

# Also check if we're in GitHub Actions or if HOOMD is not available
in_ci = (os.environ.get('CI') == 'true' or 
         os.environ.get('GITHUB_ACTIONS') == 'true')

# Try to import HOOMD to see if it's available
try:
    import hoomd as _test_hoomd
    hoomd_available = True
    del _test_hoomd
except ImportError:
    hoomd_available = False

# Need mocking if on RTD, in CI, or HOOMD not available
needs_mocking = on_rtd or in_ci or not hoomd_available

print(f"On RTD: {on_rtd}")
print(f"In CI: {in_ci}")
print(f"HOOMD available: {hoomd_available}")
print(f"Needs mocking: {needs_mocking}")

if needs_mocking:
    print("Setting up simplified plugin imports for documentation build...")
    
    # Mock only the C++ extensions (the compiled .so files)
    from unittest.mock import MagicMock
    sys.modules['hoomd.cavitymd._cavitymd'] = MagicMock()
    sys.modules['hoomd.bussi_reservoir._bussi_reservoir'] = MagicMock()
    sys.modules['cavitymd._cavitymd'] = MagicMock()
    sys.modules['bussi_reservoir._bussi_reservoir'] = MagicMock()
    print(" Mocked C++ extensions")
    
    # Try to import HOOMD (available via conda on RTD, but not in GitHub Actions)
    try:
        import hoomd
        print(" HOOMD base package imported from conda")
    except ImportError:
        # HOOMD not available - create a minimal mock
        print(" HOOMD not available - creating mock base package")
        from types import ModuleType
        
        # Create HOOMD as a package (not just a module)
        hoomd = ModuleType('hoomd')
        hoomd.__package__ = 'hoomd'
        hoomd.__path__ = ['hoomd']
        
        # Add version attribute that plugins check
        hoomd_version = ModuleType('hoomd.version')
        hoomd_version.version = '4.8.2'  # Match the required version
        hoomd.version = hoomd_version
        
        # Pre-register hoomd in sys.modules before creating submodules
        sys.modules['hoomd'] = hoomd
        sys.modules['hoomd.version'] = hoomd_version
        
        # Create bussi_reservoir submodule structure
        hoomd_bussi = ModuleType('hoomd.bussi_reservoir')
        hoomd_bussi.__package__ = 'hoomd.bussi_reservoir'
        hoomd_bussi_thermostats = ModuleType('hoomd.bussi_reservoir.thermostats')
        hoomd_bussi_thermostats.BussiReservoir = type('BussiReservoir', (), {})
        hoomd_bussi.thermostats = hoomd_bussi_thermostats
        hoomd.bussi_reservoir = hoomd_bussi
        sys.modules['hoomd.bussi_reservoir'] = hoomd_bussi
        sys.modules['hoomd.bussi_reservoir.thermostats'] = hoomd_bussi_thermostats
        
        print(" Created mock HOOMD base package with version info and bussi_reservoir")
    
    # Import our plugins directly and register them in the hoomd namespace
    try:
        # Set environment variable to bypass version check in plugin
        os.environ['RTD_ENV_NAME'] = 'docs'
        os.environ['READTHEDOCS'] = 'True'
        
        import cavitymd
        import bussi_reservoir
        
        # Register plugins in the hoomd namespace
        hoomd.cavitymd = cavitymd
        hoomd.bussi_reservoir = bussi_reservoir
        
        # Also register them in sys.modules for autosummary
        sys.modules['hoomd.cavitymd'] = cavitymd
        sys.modules['hoomd.bussi_reservoir'] = bussi_reservoir
        
        # Register all submodules and their classes for autosummary
        for module_name in ['analysis', 'forces', 'simulation', 'utils', 'variants', 'updaters']:
            if hasattr(cavitymd, module_name):
                submodule = getattr(cavitymd, module_name)
                # Register the submodule
                setattr(hoomd.cavitymd, module_name, submodule)
                sys.modules[f'hoomd.cavitymd.{module_name}'] = submodule
                
                # Also register individual classes within each submodule
                for attr_name in dir(submodule):
                    attr = getattr(submodule, attr_name)
                    # Register classes and functions (not internal attributes)
                    if not attr_name.startswith('_') and (
                        (hasattr(attr, '__call__') and hasattr(attr, '__name__')) or  # Functions
                        (hasattr(attr, '__module__') and not callable(attr)) or        # Classes
                        type(attr).__name__ in ['type', 'ABCMeta']                     # Class types
                    ):
                        # Make the class available through the hoomd namespace
                        if not hasattr(getattr(hoomd.cavitymd, module_name), attr_name):
                            setattr(getattr(hoomd.cavitymd, module_name), attr_name, attr)
        
        # Special handling for commonly used classes - make them directly accessible
        # This ensures classes like CavityForce are available as hoomd.cavitymd.forces.CavityForce
        if hasattr(cavitymd.forces, 'CavityForce'):
            hoomd.cavitymd.forces.CavityForce = cavitymd.forces.CavityForce
        if hasattr(cavitymd.simulation, 'CavityMDSimulation'):
            hoomd.cavitymd.simulation.CavityMDSimulation = cavitymd.simulation.CavityMDSimulation
        if hasattr(cavitymd.simulation, 'AdaptiveTimestepUpdater'):
            hoomd.cavitymd.simulation.AdaptiveTimestepUpdater = cavitymd.simulation.AdaptiveTimestepUpdater
        if hasattr(bussi_reservoir, 'BussiReservoir'):
            hoomd.bussi_reservoir.BussiReservoir = bussi_reservoir.BussiReservoir
            
        print(" Successfully imported and registered plugin modules")
        print(f"  CavityForce available: {hasattr(hoomd.cavitymd.forces, 'CavityForce') if hasattr(hoomd.cavitymd, 'forces') else False}")
        print(f"  BussiReservoir available: {hasattr(hoomd.bussi_reservoir, 'BussiReservoir')}")
        print(f"  Analysis module available: {hasattr(hoomd.cavitymd, 'analysis')}")
        print(f"  Simulation module available: {hasattr(hoomd.cavitymd, 'simulation')}")
        print(f"  AdaptiveTimestepUpdater available: {hasattr(hoomd.cavitymd.simulation, 'AdaptiveTimestepUpdater') if hasattr(hoomd.cavitymd, 'simulation') else False}")
        
    except Exception as e:
        print(f" Failed to import plugins: {e}")
        import traceback
        traceback.print_exc()
        
        # Create fallback empty modules to prevent autosummary errors
        from types import ModuleType
        
        # Create mock plugin modules
        mock_cavitymd = ModuleType('cavitymd')
        mock_bussi = ModuleType('bussi_reservoir')
        
        # Create mock submodules and classes
        for submodule in ['analysis', 'forces', 'simulation', 'utils', 'variants', 'updaters', 'controllers', 'data', 'experiment', 'state_manager']:
            mock_submodule = ModuleType(f'cavitymd.{submodule}')
            setattr(mock_cavitymd, submodule, mock_submodule)
            sys.modules[f'hoomd.cavitymd.{submodule}'] = mock_submodule
            
            # Add mock classes for each submodule
            if submodule == 'forces':
                mock_submodule.CavityForce = type('CavityForce', (), {})
                mock_submodule.DipoleResponseForce = type('DipoleResponseForce', (), {})
                mock_submodule.PerturbationForce = type('PerturbationForce', (), {})
            elif submodule == 'simulation':
                mock_submodule.CavityMDSimulation = type('CavityMDSimulation', (), {})
                mock_submodule.AdaptiveTimestepUpdater = type('AdaptiveTimestepUpdater', (), {})
            elif submodule == 'analysis':
                mock_submodule.Status = type('Status', (), {})
                mock_submodule.ElapsedTimeTracker = type('ElapsedTimeTracker', (), {})
                mock_submodule.TimestepFormatter = type('TimestepFormatter', (), {})
                mock_submodule.AutocorrelationTracker = type('AutocorrelationTracker', (), {})
                mock_submodule.FieldAutocorrelationTracker = type('FieldAutocorrelationTracker', (), {})
                mock_submodule.EnergyTracker = type('EnergyTracker', (), {})
                mock_submodule.DipoleAutocorrelation = type('DipoleAutocorrelation', (), {})
                mock_submodule.DipoleMomentFDRTracker = type('DipoleMomentFDRTracker', (), {})
                mock_submodule.PerformanceTracker = type('PerformanceTracker', (), {})
                mock_submodule.TemperatureTracker = type('TemperatureTracker', (), {})
                mock_submodule.CavityModeTracker = type('CavityModeTracker', (), {})
            elif submodule == 'variants':
                mock_submodule.StepVariant = type('StepVariant', (), {})
                mock_submodule.ConstantVariant = type('ConstantVariant', (), {})
            elif submodule == 'updaters':
                mock_submodule.CavityParticleDisplacer = type('CavityParticleDisplacer', (), {})
            elif submodule == 'controllers':
                mock_submodule.AdaptiveMPCController = type('AdaptiveMPCController', (), {})
                mock_submodule.DiffEqController = type('DiffEqController', (), {})
                mock_submodule.PIDControl = type('PIDControl', (), {})
                mock_submodule.SimpleSetpointController = type('SimpleSetpointController', (), {})
            elif submodule == 'data':
                mock_submodule.ObservableWriter = type('ObservableWriter', (), {})
                mock_submodule.ObservableReader = type('ObservableReader', (), {})
            elif submodule == 'experiment':
                mock_submodule.ExperimentManager = type('ExperimentManager', (), {})
                mock_submodule.ParameterSweep = type('ParameterSweep', (), {})
            elif submodule == 'state_manager':
                mock_submodule.StateManager = type('StateManager', (), {})
                mock_submodule.StateValidator = type('StateValidator', (), {})
            elif submodule == 'utils':
                mock_submodule.PhysicalConstants = type('PhysicalConstants', (), {})
                mock_submodule.unwrap_positions = lambda x, y, z: x
                mock_submodule.get_slurm_info = lambda: {}
                mock_submodule.parse_replicas = lambda x: [1]
        
        # Add top-level mock classes (these appear directly under cavitymd)
        mock_cavitymd.CompositeVariant = type('CompositeVariant', (), {})
        mock_cavitymd.FDRIntegration = type('FDRIntegration', (), {})
        mock_cavitymd.FDRTemperatureEstimator = type('FDRTemperatureEstimator', (), {})
        mock_cavitymd.FDRWorkflow = type('FDRWorkflow', (), {})
        
        # Add mock BussiReservoir
        mock_bussi.BussiReservoir = type('BussiReservoir', (), {})
        
        hoomd.cavitymd = mock_cavitymd
        hoomd.bussi_reservoir = mock_bussi
        sys.modules['hoomd.cavitymd'] = mock_cavitymd
        sys.modules['hoomd.bussi_reservoir'] = mock_bussi
        
        print(" Created fallback mock modules for documentation build")
    
    print(" Plugin setup complete")
else:
    print("Local environment with HOOMD installed - no mocking needed")
    print("Plugins will be imported directly from installed HOOMD")

# -- Path setup (redundant but kept for clarity) ------------------------------

# Add the project root to the Python path
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

# Suppress specific warnings
suppress_warnings = [
    'autodoc.duplicate_object',
    'autosummary.duplicate_object',
    'ref.duplicate',
    'toc.not_included',
    'duplicate_object',
    'py.duplicate_object',
    'duplicate',
]

# Autodoc configuration to handle duplicates better
autodoc_member_order = 'bysource'
autodoc_inherit_docstrings = True
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# Additional configuration to handle the duplicate object issue
nitpicky = False

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
    'README_DEPLOYMENT.md',
]

# The suffix(es) of source filenames.
# Let extensions auto-register their parsers
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
    '.ipynb': 'nbsphinx',
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
        "color-background-primary": "#1a1a1a",
        "color-background-secondary": "#2d2d2d",
        "color-background-hover": "#404040",
        "color-background-hover--transparent": "#40404000",
        "color-background-border": "#404040",
        
        # Foreground colors
        "color-foreground-primary": "#ffffff",
        "color-foreground-secondary": "#b0b0b0",
        "color-foreground-muted": "#808080",
        "color-foreground-border": "#666666",
        
        # API documentation colors
        "color-api-background": "#2d2d2d",
        "color-api-background-hover": "#404040",
        "color-api-overall": "#3498db",
        "color-api-name": "#ecf0f1",
        "color-api-pre-name": "#bdc3c7",
        "color-api-paren": "#95a5a6",
        "color-api-keyword": "#bb86fc",
        
        # Admonition colors
        "color-admonition-background": "transparent",
        
        # Inline code colors
        "color-inline-code-background": "#2d2d2d",
        
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
    'undoc-members': True,
    'show-inheritance': True,
    'member-order': 'bysource',
    'special-members': '__init__',
}

autodoc_typehints = 'description'
autodoc_typehints_description_target = 'documented'

# Don't fail on import errors and handle gracefully
autodoc_preserve_defaults = True

# Additional autodoc configuration for better documentation
autodoc_member_order = 'bysource'
autodoc_class_signature = 'mixed'
autodoc_type_aliases = {}

# Mock problematic imports
autodoc_mock_imports = []

# Generate autosummary files automatically
autosummary_generate = True
autosummary_imported_members = True

# -- Options for intersphinx extension ---------------------------------------

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'hoomd': ('https://hoomd-blue.readthedocs.io/en/stable/', None),
}

# -- Options for todo extension ---------------------------------------------

todo_include_todos = True

# -- Options for coverage extension ------------------------------------------

coverage_show_missing_items = True

# -- Options for Napoleon extension ------------------------------------------

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
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# -- Options for bibliographic references ------------------------------------

bibtex_bibfiles = ['references.bib']
bibtex_default_style = 'unsrt'

# -- Options for nbsphinx extension ------------------------------------------

nbsphinx_execute = 'never'  # Don't execute notebooks during build
nbsphinx_allow_errors = True
nbsphinx_kernel_name = 'python3'

# Custom CSS for better notebook rendering
nbsphinx_prolog = """
{% set docname = env.docname.replace("examples/", "") %}
{% if env.metadata[docname] %}
{% set nbpath = env.metadata[docname]["nbpath"] %}
{% else %}
{% set nbpath = env.docname + ".ipynb" %}
{% endif %}

.. raw:: html

    <div class="admonition note">
      <p class="admonition-title">Note</p>
      <p>This page was generated from a Jupyter notebook that can be 
         found in the examples directory of the project.</p>
    </div>
"""

# -- Options for myst parser ------------------------------------------------

myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]

# -- Custom setup function --------------------------------------------------

def setup(app):
    """Custom setup function for Sphinx."""
    # Add custom CSS for better API documentation styling
    app.add_css_file('custom.css')
    
    # Set up the documentation build environment
    if needs_mocking:
        print("Sphinx setup: Running in CI/RTD environment with mocking")
    else:
        print("Sphinx setup: Running locally with full HOOMD installation") 