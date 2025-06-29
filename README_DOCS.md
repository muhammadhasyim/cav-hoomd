# Cavity HOOMD Documentation System

This directory contains a comprehensive documentation website for Cavity HOOMD, built with Sphinx and designed for modern web deployment.

## 🚀 Quick Setup

### Automated Setup (Recommended)

Run the automated setup script to get everything configured:

```bash
# For GitHub Pages deployment
python setup_docs.py --deploy-method github-pages --github-repo yourusername/cavity-hoomd

# For ReadTheDocs deployment  
python setup_docs.py --deploy-method readthedocs --rtd-project cavity-hoomd

# For local development only
python setup_docs.py --deploy-method local
```

### Manual Setup

If you prefer manual setup:

```bash
# Install dependencies
pip install -r docs/requirements.txt

# Build documentation
cd docs
make html

# Serve locally
make serve
```

## 📚 Documentation Features

### Professional Website
- **Modern Design**: Clean, responsive design with dark/light mode support
- **Fast Search**: Full-text search across all documentation
- **Mobile Friendly**: Optimized for all device sizes
- **Print Support**: Clean print layouts for offline reading

### Comprehensive Content
- **API Reference**: Auto-generated from docstrings with cross-references
- **User Guides**: Step-by-step tutorials and best practices
- **Examples**: Working code examples and Jupyter notebooks
- **Theory**: Mathematical background and scientific context

### Advanced Features
- **Live Reload**: Development server with automatic rebuilds
- **Multiple Formats**: HTML, PDF, and EPUB generation
- **Code Highlighting**: Syntax highlighting for Python, bash, and more
- **Math Support**: LaTeX math rendering with MathJax
- **Citations**: BibTeX bibliography support
- **Copy Buttons**: One-click code copying

## 🏗️ Project Structure

```
docs/
├── conf.py                 # Sphinx configuration
├── index.rst              # Homepage
├── quickstart.rst         # Quick start guide
├── requirements.txt       # Python dependencies
├── Makefile               # Build commands
│
├── api/                   # API reference
│   ├── index.rst
│   ├── forces.rst
│   ├── simulation.rst
│   ├── analysis.rst
│   └── utils.rst
│
├── user_guide/            # User guides  
│   ├── index.rst
│   ├── basic_usage.rst
│   ├── advanced_features.rst
│   ├── parameter_sweeps.rst
│   ├── analysis.rst
│   └── performance.rst
│
├── examples/              # Examples and tutorials
│   ├── index.rst
│   ├── basic_simulation.ipynb
│   ├── energy_tracking.ipynb
│   └── parameter_sweep.ipynb
│
├── _static/              # Static files
│   ├── custom.css        # Custom styling
│   ├── logo.png          # Project logo
│   └── favicon.ico       # Browser icon
│
└── _build/               # Generated files (auto-created)
    └── html/             # Built website
```

## 🛠️ Build Commands

The documentation uses a Makefile with convenient commands:

```bash
cd docs

# Development
make live          # Live reload server (auto-rebuilds)
make html          # Build HTML documentation
make serve         # Serve built documentation

# Quality control  
make check         # Check for warnings and errors
make linkcheck     # Verify external links

# Multiple formats
make pdf           # Build PDF (requires LaTeX)
make epub          # Build EPUB

# Maintenance
make clean-all     # Clean all build files
make dev           # Development workflow (clean + build)
make prod          # Production workflow (clean + check + build)
```

## 🌐 Deployment Options

### 1. GitHub Pages (Recommended)

Automatic deployment with GitHub Actions:

```bash
# Setup (run once)
python setup_docs.py --deploy-method github-pages --github-repo username/repo

# Deploy
git add .
git commit -m "Add documentation"
git push origin main

# Enable GitHub Pages in repository settings
# Set source to "GitHub Actions"
# Visit: https://username.github.io/repo/
```

**Features:**
- ✅ Free hosting
- ✅ Automatic builds on push
- ✅ Custom domains supported
- ✅ Built-in SSL

### 2. ReadTheDocs

Professional documentation hosting:

```bash
# Setup
python setup_docs.py --deploy-method readthedocs --rtd-project cavity-hoomd

# Deploy
git push origin main

# Import project at readthedocs.org
# Visit: https://cavity-hoomd.readthedocs.io/
```

**Features:**
- ✅ Free for open source
- ✅ PDF/EPUB downloads
- ✅ Versioning support
- ✅ Analytics included
- ✅ Search optimization

### 3. Netlify

Modern web hosting with edge CDN:

```bash
# Setup
python setup_docs.py --deploy-method netlify

# Deploy
git push origin main

# Connect repository at netlify.com
# Automatic deployment on commits
```

**Features:**
- ✅ Fast global CDN
- ✅ Branch previews
- ✅ Form handling
- ✅ Serverless functions

### 4. Local Development

For development and testing:

```bash
cd docs
make live     # Start development server
# Visit: http://localhost:8000
```

## ✏️ Content Creation

### Adding New Pages

1. **Create RST file**: Add new `.rst` file in appropriate directory
2. **Update TOC**: Add to `toctree` directive in parent `index.rst`  
3. **Build**: Run `make html` to generate

Example:
```rst
My New Page
===========

Content goes here with **bold** and *italic* text.

Code Examples
-------------

.. code-block:: python

   # Python code example
   import hoomd
   print("Hello, World!")

Math Equations
--------------

.. math::

   E = mc^2
```

### Adding API Documentation

API docs are auto-generated from docstrings:

```python
class MyClass:
    """One-line description.
    
    Longer description with more details about the class
    functionality and usage patterns.
    
    Parameters
    ----------
    param1 : str
        Description of param1
    param2 : int, optional
        Description of param2 (default: 42)
        
    Examples
    --------
    >>> obj = MyClass("hello", 123)
    >>> obj.do_something()
    """
    
    def do_something(self):
        """Do something useful.
        
        Returns
        -------
        bool
            True if successful
        """
        return True
```

### Adding Examples

Create Jupyter notebooks in `examples/`:

```bash
# Create new notebook
jupyter notebook examples/my_example.ipynb

# Add to examples/index.rst toctree:
# my_example.ipynb
```

## 🎨 Customization

### Styling

Edit `docs/_static/custom.css` for visual customization:

```css
/* Brand colors */
:root {
    --color-brand-primary: #336790;
    --color-brand-content: #336790;
}

/* Custom styles */
.my-custom-class {
    background-color: #f8f9fa;
    padding: 1rem;
    border-radius: 0.5rem;
}
```

### Configuration

Edit `docs/conf.py` for Sphinx settings:

```python
# Project information
project = 'Your Project Name'
author = 'Your Name'

# Theme options
html_theme_options = {
    "source_repository": "https://github.com/username/repo",
    "source_branch": "main",
}

# Extensions
extensions = [
    'sphinx.ext.autodoc',
    'sphinx_design',
    # Add more extensions
]
```

### Logo and Branding

Replace files in `docs/_static/`:
- `logo.png` - Main logo (200x100px recommended)
- `favicon.ico` - Browser icon (32x32px)

## 📊 Analytics and Monitoring

### GitHub Pages
- Use Google Analytics by adding tracking code to `_templates/layout.html`
- Monitor traffic in GitHub Insights

### ReadTheDocs
- Built-in analytics available in dashboard
- Traffic stats and search analytics included

### Netlify
- Built-in analytics in Netlify dashboard
- Real-time visitor monitoring

## 🔧 Troubleshooting

### Common Issues

**Build Fails with Import Errors**
```bash
# Solution: Install your package in development mode
pip install -e .
cd docs && make html
```

**Math Not Rendering**
```bash
# Check MathJax configuration in conf.py
mathjax3_config = {
    'tex': {'inlineMath': [['$', '$'], ['\\(', '\\)']]},
}
```

**Links Broken**
```bash
# Check for broken external links
cd docs && make linkcheck
```

**Slow Builds**
```bash
# Use parallel builds
make html SPHINXOPTS="-j auto"
```

### Getting Help

1. **Check Build Logs**: Look at build output for specific errors
2. **Validate RST**: Use online RST validators for syntax checking
3. **Sphinx Documentation**: Comprehensive guide at sphinx-doc.org
4. **GitHub Issues**: Report problems or ask questions

## 📈 Performance Optimization

### Build Speed
- Use `sphinx-autobuild` for development
- Enable parallel builds with `-j auto`
- Exclude large files from builds

### Page Load Speed
- Optimize images (use WebP format)
- Minimize custom CSS
- Enable CDN for static assets

### SEO Optimization
- Add meta descriptions to pages
- Use semantic HTML structure
- Include sitemap.xml (auto-generated)

## 🔄 Maintenance

### Regular Tasks

**Weekly**
- Check for broken external links: `make linkcheck`
- Update dependencies: `pip install -U -r requirements.txt`
- Review build warnings and fix issues

**Monthly**
- Update Sphinx and theme versions
- Review and update content for accuracy
- Check analytics for popular/problematic pages

**Version Releases**
- Tag documentation versions
- Update changelog
- Archive old documentation if needed

### Automation

Set up automated checks:

```yaml
# .github/workflows/docs-check.yml
name: Documentation Check
on:
  pull_request:
    paths: ['docs/**']

jobs:
  check:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Check documentation
        run: |
          cd docs
          make check
          make linkcheck
```

## 🤝 Contributing

### For Writers
1. Fork the repository
2. Create a feature branch
3. Add/edit documentation
4. Test builds locally
5. Submit pull request

### For Developers
1. Add docstrings to new code
2. Update API docs when changing interfaces
3. Add examples for new features
4. Test documentation builds

### Style Guide
- Use clear, concise language
- Include code examples
- Add cross-references with `:doc:`, `:class:`, etc.
- Follow existing formatting patterns

---

## 📝 License

Documentation is licensed under the same terms as Cavity HOOMD (BSD 3-Clause License).

## 🙏 Acknowledgments

Built with:
- [Sphinx](https://www.sphinx-doc.org/) - Documentation generator
- [Furo](https://pradyunsg.me/furo/) - Modern theme
- [ReadTheDocs](https://readthedocs.org/) - Hosting platform
- [GitHub Pages](https://pages.github.com/) - Free hosting 