# Cavity HOOMD Documentation

Simple, focused documentation for the Cavity HOOMD molecular dynamics plugin.

## 🚀 Quick Setup

Build the documentation:

```bash
cd docs
pip install -r requirements.txt
make html
```

The documentation will be in `docs/_build/html/`.

## 📚 Documentation Structure

The documentation is now streamlined and focused on the single main example:

```
docs/
├── index.rst              # Main documentation page
├── quickstart.rst         # Quick start guide
├── installation.rst       # Installation instructions
├── theory.rst             # Simplified theory background
├── license.rst            # License information
├── conf.py               # Sphinx configuration
├── requirements.txt      # Documentation dependencies
└── Makefile             # Build commands
```

## 🎯 Focus

The documentation centers around the `examples/05_advanced_run.py` script, which provides:

- Complete command-line interface for cavity MD simulations
- Built-in parameter sweeps and replica management
- Comprehensive output and analysis

## 🛠️ Build Commands

```bash
cd docs

# Build HTML documentation
make html

# Serve locally
make serve

# Clean build
make clean
```

## 🌐 Deployment

### GitHub Pages

1. Enable GitHub Pages in repository settings
2. Set source to "Deploy from a branch" 
3. Select `gh-pages` branch
4. Push changes to trigger automatic deployment

### ReadTheDocs

1. Import project at readthedocs.org
2. Link to your GitHub repository
3. Documentation builds automatically on commits

## 📋 Content Guidelines

The simplified documentation focuses on:

1. **Getting Started** - Quick installation and first simulation
2. **Essential Usage** - Core command-line options and examples  
3. **Basic Theory** - Accessible physics background
4. **Reference** - Installation and licensing information

This streamlined approach makes it easier for users to quickly understand and use the software without getting overwhelmed by excessive detail. 