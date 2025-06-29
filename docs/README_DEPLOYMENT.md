# Documentation Deployment Guide

This guide explains how to deploy the Cavity HOOMD documentation to GitHub Pages with automatic builds and HOOMD-blue styling.

## Overview

The documentation system is configured with:
- **Sphinx** for documentation generation
- **Furo theme** with HOOMD-blue styling
- **GitHub Actions** for automatic deployment
- **GitHub Pages** for hosting

## Quick Setup

### 1. Enable GitHub Pages

1. Go to your repository: https://github.com/muhammadhasyim/cav-hoomd
2. Click **Settings** → **Pages**
3. Under **Source**, select **"GitHub Actions"**
4. Save the settings

### 2. Push Changes

The documentation will automatically build and deploy when you:
- Push to the `main` branch
- Modify files in `docs/`, `src/`, `examples/`, or documentation files

### 3. Access Your Documentation

After the first successful deployment, your documentation will be available at:
```
https://muhammadhasyim.github.io/cav-hoomd/
```

## Manual Deployment

You can also trigger a manual deployment:

1. Go to **Actions** tab in your GitHub repository
2. Click **Documentation** workflow
3. Click **Run workflow**
4. Check **"Deploy to GitHub Pages"**
5. Click **Run workflow**

## Local Development

### Prerequisites

Install the documentation dependencies:
```bash
cd docs
pip install -r requirements.txt
```

### Build Locally

```bash
cd docs
make html
```

The built documentation will be in `docs/_build/html/`.

### Live Preview

For live reload during development:
```bash
cd docs
make live
```

This will start a local server at `http://localhost:8000` with automatic reload.

## Styling and Theme

The documentation uses the **Furo theme** configured to match **HOOMD-blue**'s appearance:

### Color Scheme
- **Primary brand color**: `#2980b9` (HOOMD-blue inspired)
- **Light mode**: Clean whites and grays
- **Dark mode**: Dark backgrounds with blue accents
- **Code highlighting**: Syntax-aware with appropriate contrast

### Typography
- **Font**: Inter (fallback to system fonts)
- **Headings**: Weighted hierarchy with brand color accents
- **Code**: Monospace with syntax highlighting

### Features
- Light/dark mode toggle
- Responsive design
- Mobile-friendly navigation
- Copy buttons for code blocks
- Search functionality
- Edit page buttons linked to GitHub

## Customization

### Theme Colors

Edit `docs/conf.py` to modify the color scheme:
```python
html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": "#your-color",
        # ... other colors
    },
    "dark_css_variables": {
        "color-brand-primary": "#your-dark-color",
        # ... other colors
    },
}
```

### Custom CSS

Add custom styles to `docs/_static/custom.css`. The file includes:
- HOOMD-blue inspired styling
- API documentation improvements
- Navigation enhancements
- Responsive design improvements

### Logo and Favicon

Place your logo and favicon in `docs/_static/`:
- `logo.png` - Main logo (recommended: 200px wide)
- `favicon.ico` - Browser favicon (16x16 or 32x32 pixels)

## Workflow Configuration

The GitHub Actions workflow (`.github/workflows/docs.yaml`) includes:

### Triggers
- Push to `main`/`master` branches
- Pull requests (build only)
- Manual dispatch

### Jobs
1. **Build**: Installs dependencies, builds documentation
2. **Deploy**: Uploads to GitHub Pages (main branch only)
3. **Link Check**: Validates links on PRs

### Environment Variables
- `PYTHONPATH`: Includes project root and source directories
- Build options: Warnings as errors, keep going, nitpicky mode

## Troubleshooting

### Common Issues

1. **Build Failures**
   - Check the Actions tab for error logs
   - Ensure all dependencies are in `requirements.txt`
   - Verify Python imports work correctly

2. **Missing Files**
   - Ensure `docs/_static/` and `docs/_templates/` exist
   - Check that referenced files are committed to git

3. **Styling Issues**
   - Clear browser cache after CSS changes
   - Check `custom.css` for syntax errors
   - Verify CSS variable names match Furo's expectations

4. **GitHub Pages Not Working**
   - Verify Pages is enabled in repository settings
   - Check that the workflow has proper permissions
   - Ensure the deployment job completed successfully

### Debug Mode

For verbose build output, modify the Makefile or run:
```bash
cd docs
sphinx-build -W --keep-going -n -v -b html . _build/html
```

## Advanced Features

### API Documentation

The system automatically generates API documentation from docstrings using:
- `sphinx.ext.autodoc`
- `sphinx.ext.autosummary`
- `sphinx_autodoc_typehints`

### Jupyter Notebooks

Include Jupyter notebooks with `nbsphinx`:
1. Place `.ipynb` files in appropriate directories
2. Reference them in `.rst` files with `.. toctree::`

### Mathematical Expressions

Use MathJax for equations:
```rst
.. math::
   
   E = mc^2
```

Or inline: :math:`F = ma`

### Bibliography

Add citations using `sphinxcontrib.bibtex`:
1. Create `references.bib` in the docs directory
2. Use `:cite:` roles to reference entries

## Maintenance

### Regular Updates

1. **Dependencies**: Update `requirements.txt` monthly
2. **Theme**: Check for Furo theme updates quarterly
3. **Actions**: Update GitHub Actions versions annually

### Monitoring

- Monitor build times and success rates
- Check for broken links periodically  
- Review documentation coverage

### Backup

The documentation source is in git, but consider:
- Backing up built HTML periodically
- Documenting any custom extensions or modifications
- Maintaining deployment configuration as code

## Support

For issues with:
- **Sphinx**: Check [Sphinx documentation](https://www.sphinx-doc.org/)
- **Furo theme**: Check [Furo documentation](https://pradyunsg.me/furo/)
- **GitHub Actions**: Check [GitHub Actions documentation](https://docs.github.com/en/actions)
- **This setup**: Open an issue in this repository 