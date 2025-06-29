#!/usr/bin/env python3
"""
Automated setup script for Cavity HOOMD documentation website.

This script sets up a complete documentation system including:
- Sphinx configuration
- GitHub Actions for automated builds
- ReadTheDocs configuration
- Local development environment

Usage:
    python setup_docs.py [--deploy-method METHOD] [--github-repo REPO] [--rtd-project PROJECT]

Deploy methods:
    - github-pages: Deploy to GitHub Pages
    - readthedocs: Deploy to ReadTheDocs
    - netlify: Deploy to Netlify
    - local: Local development only
"""

import os
import sys
import argparse
import shutil
from pathlib import Path
import subprocess


class DocumentationSetup:
    """Setup and configure documentation system for Cavity HOOMD."""
    
    def __init__(self, deploy_method='github-pages', github_repo=None, rtd_project=None):
        self.deploy_method = deploy_method
        self.github_repo = github_repo or 'yourusername/cavity-hoomd'
        self.rtd_project = rtd_project or 'cavity-hoomd'
        self.docs_dir = Path('docs')
        self.project_root = Path('.')
        
    def setup_all(self):
        """Run complete documentation setup."""
        print("🚀 Setting up Cavity HOOMD Documentation Website")
        print("=" * 60)
        
        # Check prerequisites
        self.check_prerequisites()
        
        # Create directory structure
        self.create_directory_structure()
        
        # Setup deployment configuration
        if self.deploy_method == 'github-pages':
            self.setup_github_pages()
        elif self.deploy_method == 'readthedocs':
            self.setup_readthedocs()
        elif self.deploy_method == 'netlify':
            self.setup_netlify()
        
        # Create additional files
        self.create_additional_files()
        
        # Install dependencies
        self.install_dependencies()
        
        # Generate initial build
        self.build_documentation()
        
        print("\n✅ Documentation setup complete!")
        self.print_next_steps()
    
    def check_prerequisites(self):
        """Check that required tools are available."""
        print("📋 Checking prerequisites...")
        
        # Check Python
        if sys.version_info < (3, 8):
            raise RuntimeError("Python 3.8+ required")
        
        # Check if git is available
        try:
            subprocess.run(['git', '--version'], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("⚠️  Git not found - GitHub deployment will not work")
        
        # Check if we're in a git repository
        if not (self.project_root / '.git').exists():
            print("⚠️  Not in a git repository - some features may not work")
        
        print("✅ Prerequisites checked")
    
    def create_directory_structure(self):
        """Create the basic directory structure."""
        print("📁 Creating directory structure...")
        
        # Create main directories
        directories = [
            'docs/_static',
            'docs/_templates', 
            'docs/api',
            'docs/user_guide',
            'docs/examples',
            'docs/development',
            'docs/images',
        ]
        
        for directory in directories:
            Path(directory).mkdir(parents=True, exist_ok=True)
        
        print("✅ Directory structure created")
    
    def setup_github_pages(self):
        """Setup GitHub Pages deployment."""
        print("🐙 Setting up GitHub Pages deployment...")
        
        # Create GitHub Actions workflow
        github_dir = Path('.github/workflows')
        github_dir.mkdir(parents=True, exist_ok=True)
        
        workflow_content = f"""name: Build and Deploy Documentation

on:
  push:
    branches: [ main, master ]
  pull_request:
    branches: [ main, master ]

permissions:
  contents: read
  pages: write
  id-token: write

concurrency:
  group: "pages"
  cancel-in-progress: false

jobs:
  build:
    runs-on: ubuntu-latest
    
    steps:
    - name: Checkout
      uses: actions/checkout@v4
      
    - name: Setup Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.10'
        
    - name: Cache pip packages
      uses: actions/cache@v3
      with:
        path: ~/.cache/pip
        key: ${{{{ runner.os }}}}-pip-${{{{ hashFiles('docs/requirements.txt') }}}}
        
    - name: Install dependencies
      run: |
        pip install -r docs/requirements.txt
        
    - name: Build documentation
      run: |
        cd docs
        make html
        
    - name: Setup Pages
      if: github.ref == 'refs/heads/main' && github.event_name == 'push'
      uses: actions/configure-pages@v3
      
    - name: Upload artifact
      if: github.ref == 'refs/heads/main' && github.event_name == 'push'
      uses: actions/upload-pages-artifact@v2
      with:
        path: docs/_build/html

  deploy:
    if: github.ref == 'refs/heads/main' && github.event_name == 'push'
    environment:
      name: github-pages
      url: ${{{{ steps.deployment.outputs.page_url }}}}
    runs-on: ubuntu-latest
    needs: build
    steps:
      - name: Deploy to GitHub Pages
        id: deployment
        uses: actions/deploy-pages@v2
"""
        
        with open(github_dir / 'docs.yml', 'w') as f:
            f.write(workflow_content)
        
        print(f"✅ GitHub Pages workflow created")
        print(f"   Repository: {self.github_repo}")
        print(f"   Documentation will be available at: https://{self.github_repo.split('/')[0]}.github.io/{self.github_repo.split('/')[1]}/")
    
    def setup_readthedocs(self):
        """Setup ReadTheDocs deployment."""
        print("📖 Setting up ReadTheDocs deployment...")
        
        # Create .readthedocs.yaml
        rtd_config = f"""# .readthedocs.yaml
# Read the Docs configuration file
# See https://docs.readthedocs.io/en/stable/config-file/v2.html for details

version: 2

build:
  os: ubuntu-22.04
  tools:
    python: "3.10"

sphinx:
  configuration: docs/conf.py
  builder: html
  fail_on_warning: true

formats:
  - pdf
  - epub

python:
  install:
    - requirements: docs/requirements.txt
    - method: pip
      path: .
      extra_requirements:
        - docs

submodules:
  include: all
  recursive: true
"""
        
        with open('.readthedocs.yaml', 'w') as f:
            f.write(rtd_config)
        
        print(f"✅ ReadTheDocs configuration created")
        print(f"   Project: {self.rtd_project}")
        print(f"   Documentation will be available at: https://{self.rtd_project}.readthedocs.io/")
    
    def setup_netlify(self):
        """Setup Netlify deployment."""
        print("🌐 Setting up Netlify deployment...")
        
        # Create netlify.toml
        netlify_config = """[build]
  base = "docs/"
  publish = "_build/html/"
  command = "make html"

[build.environment]
  PYTHON_VERSION = "3.10"

[[headers]]
  for = "/*"
  [headers.values]
    X-Frame-Options = "DENY"
    X-XSS-Protection = "1; mode=block"
    X-Content-Type-Options = "nosniff"

[[redirects]]
  from = "/api/modules"
  to = "/api/index.html"
  status = 301
"""
        
        with open('netlify.toml', 'w') as f:
            f.write(netlify_config)
        
        print("✅ Netlify configuration created")
        print("   Connect your repository to Netlify for automatic deployment")
    
    def create_additional_files(self):
        """Create additional supporting files."""
        print("📝 Creating additional files...")
        
        # Create .gitignore for docs
        gitignore_content = """# Sphinx build files
_build/
_static/
_templates/

# API autosummary
api/_autosummary/

# Jupyter notebook checkpoints
**/.ipynb_checkpoints/

# OS files
.DS_Store
Thumbs.db

# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
env/
venv/
"""
        
        with open('docs/.gitignore', 'w') as f:
            f.write(gitignore_content)
        
        # Create a simple logo placeholder
        logo_placeholder = """<!-- Logo placeholder -->
<svg width="200" height="100" xmlns="http://www.w3.org/2000/svg">
  <rect width="200" height="100" fill="#336790"/>
  <text x="100" y="55" font-family="Arial, sans-serif" font-size="20" 
        fill="white" text-anchor="middle" font-weight="bold">
    Cavity HOOMD
  </text>
</svg>"""
        
        with open('docs/_static/logo.svg', 'w') as f:
            f.write(logo_placeholder)
        
        # Create bibliography file
        bibtex_content = """@article{cavity_hoomd_2025,
  title={Cavity HOOMD: A Framework for Cavity-Coupled Molecular Dynamics},
  author={Development Team},
  journal={Journal of Computational Chemistry},
  year={2025},
  volume={XX},
  pages={XXX-XXX},
  doi={10.1002/jcc.XXXXX}
}

@article{hoomd_blue_2008,
  title={HOOMD-blue: A Python Package for High-Performance Molecular Dynamics and Hard Particle Monte Carlo Simulations},
  author={Anderson, Joshua A and Glaser, Jens and Glotzer, Sharon C},
  journal={Computational Materials Science},
  volume={173},
  pages={109363},
  year={2020},
  publisher={Elsevier}
}"""
        
        with open('docs/references.bib', 'w') as f:
            f.write(bibtex_content)
        
        print("✅ Additional files created")
    
    def install_dependencies(self):
        """Install documentation dependencies."""
        print("📦 Installing documentation dependencies...")
        
        try:
            subprocess.run([
                sys.executable, '-m', 'pip', 'install', 
                '-r', 'docs/requirements.txt'
            ], check=True)
            print("✅ Dependencies installed")
        except subprocess.CalledProcessError:
            print("⚠️  Failed to install dependencies automatically")
            print("   Please run: pip install -r docs/requirements.txt")
    
    def build_documentation(self):
        """Build the documentation."""
        print("🔨 Building initial documentation...")
        
        try:
            # Change to docs directory
            os.chdir('docs')
            
            # Run make html
            subprocess.run(['make', 'html'], check=True)
            
            print("✅ Documentation built successfully")
            print(f"   Open docs/_build/html/index.html to view")
            
        except subprocess.CalledProcessError:
            print("⚠️  Failed to build documentation automatically")
            print("   Please run: cd docs && make html")
        except FileNotFoundError:
            print("⚠️  Make not found, trying direct sphinx-build")
            try:
                subprocess.run([
                    'sphinx-build', '-b', 'html', '.', '_build/html'
                ], check=True)
                print("✅ Documentation built with sphinx-build")
            except subprocess.CalledProcessError:
                print("⚠️  Failed to build with sphinx-build")
        finally:
            # Return to original directory
            os.chdir('..')
    
    def print_next_steps(self):
        """Print next steps for the user."""
        print("\n🎯 Next Steps:")
        print("=" * 30)
        
        if self.deploy_method == 'github-pages':
            print("1. Push your code to GitHub")
            print("2. Enable GitHub Pages in repository settings")
            print("3. Set source to 'GitHub Actions'")
            print(f"4. Visit https://{self.github_repo.split('/')[0]}.github.io/{self.github_repo.split('/')[1]}/")
            
        elif self.deploy_method == 'readthedocs':
            print("1. Push your code to GitHub/GitLab/Bitbucket")
            print("2. Import project at https://readthedocs.org/")
            print(f"3. Visit https://{self.rtd_project}.readthedocs.io/")
            
        elif self.deploy_method == 'netlify':
            print("1. Push your code to GitHub/GitLab/Bitbucket")
            print("2. Connect repository at https://netlify.com/")
            print("3. Netlify will auto-deploy on commits")
        
        print("\n📚 Local Development:")
        print("  cd docs")
        print("  make live      # Live reload server")
        print("  make html      # Build HTML")
        print("  make serve     # Serve built docs")
        
        print("\n🛠️  Customization:")
        print("  - Edit docs/conf.py for Sphinx settings")
        print("  - Modify docs/_static/custom.css for styling")
        print("  - Update docs/index.rst for homepage content")
        print("  - Add content to docs/user_guide/, docs/examples/, etc.")
        
        print("\n📖 Documentation Structure:")
        print("  docs/")
        print("  ├── index.rst          # Homepage")
        print("  ├── quickstart.rst     # Quick start guide")
        print("  ├── api/               # API reference") 
        print("  ├── user_guide/        # User guides")
        print("  ├── examples/          # Examples and tutorials")
        print("  └── _static/           # Static files (CSS, images)")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Setup Cavity HOOMD documentation website",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        '--deploy-method', 
        choices=['github-pages', 'readthedocs', 'netlify', 'local'],
        default='github-pages',
        help='Deployment method (default: github-pages)'
    )
    
    parser.add_argument(
        '--github-repo',
        help='GitHub repository (format: username/repo-name)'
    )
    
    parser.add_argument(
        '--rtd-project',
        help='ReadTheDocs project name'
    )
    
    args = parser.parse_args()
    
    # Create and run setup
    setup = DocumentationSetup(
        deploy_method=args.deploy_method,
        github_repo=args.github_repo,
        rtd_project=args.rtd_project
    )
    
    try:
        setup.setup_all()
    except KeyboardInterrupt:
        print("\n❌ Setup interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Setup failed: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main() 