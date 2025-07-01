# Makefile for Sphinx documentation

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD  ?= sphinx-build
SOURCEDIR    = .
BUILDDIR     = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Custom targets for Cavity HOOMD documentation

.PHONY: clean-all install-deps live serve check linkcheck

# Clean everything including autosummary files
clean-all:
	rm -rf $(BUILDDIR)/*
	rm -rf api/_autosummary/
	rm -rf examples/.ipynb_checkpoints/

# Install documentation dependencies
install-deps:
	pip install -r requirements.txt

# Build HTML documentation (single-threaded to avoid search index issues)
html:
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O) -j 1

# Build HTML documentation with warnings as errors
html-fast:
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O) -W --keep-going -j 1

# Build HTML documentation without treating warnings as errors
html-tolerant:
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O) --keep-going -j 1

# Build with live reload for development
live:
	sphinx-autobuild "$(SOURCEDIR)" "$(BUILDDIR)/html" --host 0.0.0.0 --port 8000 --open-browser

# Serve built documentation
serve:
	cd $(BUILDDIR)/html && python -m http.server 8000

# Check documentation for issues
check:
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O) -W -n --keep-going -j 1

# Check external links
linkcheck:
	@$(SPHINXBUILD) -b linkcheck "$(SOURCEDIR)" "$(BUILDDIR)/linkcheck" $(SPHINXOPTS) $(O)

# Build PDF documentation (requires LaTeX)
pdf:
	@$(SPHINXBUILD) -b latex "$(SOURCEDIR)" "$(BUILDDIR)/latex" $(SPHINXOPTS) $(O)
	@make -C $(BUILDDIR)/latex all-pdf
	@echo "PDF built in $(BUILDDIR)/latex/"

# Build EPUB documentation
epub:
	@$(SPHINXBUILD) -b epub "$(SOURCEDIR)" "$(BUILDDIR)/epub" $(SPHINXOPTS) $(O)
	@echo "EPUB built in $(BUILDDIR)/epub/"

# Generate API documentation
api:
	sphinx-apidoc -f -o api/ ../src/hoomd/cavitymd/ --separate

# Full build with all formats
all: clean-all html pdf epub
	@echo "All documentation formats built successfully!"

# Development workflow
dev: clean-all install-deps html-fast
	@echo "Development build complete!"

# Production build with checks
prod: clean-all install-deps check linkcheck html pdf
	@echo "Production build complete!"

# Show build info
info:
	@echo "Cavity HOOMD Documentation Build System"
	@echo "======================================="
	@echo "Source directory: $(SOURCEDIR)"
	@echo "Build directory:  $(BUILDDIR)"
	@echo "Sphinx build:     $(SPHINXBUILD)"
	@echo ""
	@echo "Available targets:"
	@echo "  html       - Build HTML documentation"
	@echo "  html-fast  - Fast HTML build with warnings as errors"
	@echo "  live       - Live reload development server"
	@echo "  serve      - Serve built documentation"
	@echo "  check      - Check documentation for issues"
	@echo "  linkcheck  - Check external links"
	@echo "  pdf        - Build PDF documentation"
	@echo "  epub       - Build EPUB documentation"
	@echo "  api        - Generate API documentation"
	@echo "  clean-all  - Clean all build files"
	@echo "  dev        - Development build workflow"
	@echo "  prod       - Production build workflow"
	@echo "  all        - Build all formats" 