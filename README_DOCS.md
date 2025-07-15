# Cavity HOOMD Documentation

Comprehensive, up-to-date documentation for the advanced Cavity HOOMD molecular dynamics plugin with time-varying coupling capabilities.

## 🚀 Quick Setup

Build the documentation:

```bash
cd docs
pip install -r requirements.txt
make html
```

The documentation will be in `docs/_build/html/`.

## 📚 Documentation Structure

The documentation has been completely updated to reflect all new features and capabilities:

```
docs/
├── index.rst              # Main documentation page with time-varying coupling
├── quickstart.rst         # Comprehensive quick start guide
├── installation.rst       # Updated installation with simplified architecture
├── theory.rst             # Extended theory with time-varying coupling
├── license.rst            # License information
├── conf.py               # Sphinx configuration
├── requirements.txt      # Documentation dependencies
├── api/
│   └── index.rst         # Complete API reference with all new components
└── Makefile             # Build commands
```

## 🎯 Documentation Updates

### Major Feature Additions

**Time-Varying Coupling**
- Complete mathematical formulation for dynamic coupling control
- Step-function coupling switching theory and implementation
- Energy conservation during coupling transitions
- Finite-q cavity particle displacement mechanics

**GPU Acceleration**
- Optimized CUDA kernel documentation
- Performance comparisons and benchmarks
- GPU-specific installation and usage instructions

**Advanced Analysis Framework**
- Energy conservation testing and validation
- Correlation function analysis (F(k,t))
- Cavity mode tracking and monitoring
- Performance profiling and optimization

**Simulation Framework**
- CavityMDSimulation class with comprehensive features
- Smart cavity particle handling and validation
- Parameter validation following scientific best practices
- SLURM integration for high-performance computing

**Adaptive Timestep Control**
- Automatic timestep optimization algorithms
- Stability and performance balancing
- Integration with time-varying coupling systems

### Documentation Sections Updated

**1. Main Page (index.rst)**
- Added time-varying coupling introduction
- Comprehensive feature overview
- Updated quick start examples
- Enhanced mathematical formulation

**2. Quick Start Guide (quickstart.rst)**
- Time-varying coupling examples
- Energy conservation testing procedures
- GPU acceleration usage
- Advanced analysis options
- Jupyter notebook integration
- Comprehensive parameter reference

**3. Installation Guide (installation.rst)**
- Simplified architecture overview
- Automated installation procedures
- GPU support configuration
- HPC environment setup
- Troubleshooting guide
- Development installation

**4. Theory Section (theory.rst)**
- Time-varying coupling mathematical framework
- Step-function coupling theory
- Energy conservation in dynamic systems
- Finite-q cavity particle displacement
- Dissipation in time-varying systems
- Physical interpretation of dynamic coupling

**5. API Reference (api/index.rst)**
- Complete component documentation
- Time-varying coupling variants
- Simulation framework classes
- Analysis and tracking components
- Validation and setup utilities
- Migration guide from previous versions

## 🛠️ Build Commands

```bash
cd docs

# Build HTML documentation
make html

# Serve locally for development
make serve

# Clean build artifacts
make clean

# Build and serve (development workflow)
make html && make serve
```

## ✨ Key Features Documented

### Core Capabilities
- **Time-varying coupling**: Step-function coupling switching for dynamic simulations
- **Multiple thermostats**: Bussi and Langevin thermostats with reservoir energy tracking
- **GPU acceleration**: Optimized CUDA kernels for high-performance simulations
- **Flexible coupling**: Support for both q=0 and finite-q cavity modes
- **HOOMD integration**: Native HOOMD-blue plugins with full ecosystem compatibility

### Advanced Analysis
- **Energy tracking**: Comprehensive energy component monitoring
- **Correlation functions**: F(k,t) density correlation analysis
- **Cavity mode monitoring**: Real-time cavity particle tracking
- **Performance tracking**: Detailed performance metrics and benchmarking
- **Adaptive timestep**: Automatic timestep optimization for stability

### Simulation Framework
- **Smart particle handling**: Automatic cavity particle detection and validation
- **Replica management**: Built-in support for multiple independent simulations
- **SLURM integration**: Native high-performance computing support
- **Parameter validation**: Scientific accuracy and reproducibility validation
- **Modular architecture**: Clean separation of components

### Time-Varying Experiments
- **Coupling switching**: Instantaneous coupling activation at specified times
- **Dissipation ramping**: Time-dependent cavity damping
- **Finite-q displacement**: Automatic cavity particle repositioning
- **Energy conservation**: Rigorous energy conservation testing

## 🌐 Deployment

### GitHub Pages

The documentation is automatically deployed to GitHub Pages on commits to the main branch.

### ReadTheDocs

The documentation is configured for ReadTheDocs deployment with:
- Conda environment support
- Real docstring preservation
- Automatic API documentation generation

## 📋 Content Guidelines

The updated documentation follows these principles:

1. **Scientific Accuracy** - All mathematical formulations are rigorous and correct
2. **Comprehensive Coverage** - Every new feature is fully documented with examples
3. **Practical Usage** - Real-world examples and use cases for all features
4. **Performance Focus** - Optimization guidance and performance considerations
5. **Accessibility** - Clear explanations for users at all levels

## 🔧 Architecture Improvements

The documentation now reflects the simplified plugin architecture:

**Simplified Structure**: Direct imports from `hoomd.cavitymd` and `hoomd.bussi_reservoir`
**Streamlined Documentation**: Uses conda-installed HOOMD with real docstrings
**Better Integration**: Seamless ReadTheDocs and GitHub Pages deployment
**Enhanced API**: Complete documentation of all components and their interactions

## 🧪 Testing and Validation

The documentation includes comprehensive testing guidance:

- Energy conservation validation procedures
- Time-varying coupling verification methods
- GPU vs CPU performance comparisons
- Parameter validation and error handling
- Troubleshooting common issues

## 📊 Examples and Use Cases

All documentation sections now include:

- **Basic usage examples** for quick start
- **Advanced coupling switching scenarios**
- **Energy conservation testing procedures**
- **GPU acceleration examples**
- **HPC deployment guidance**
- **Analysis and visualization workflows**

## 🔄 Migration Support

The documentation provides clear migration paths:

- From previous plugin versions
- Import statement updates
- New feature adoption guidance
- Performance optimization tips
- Best practice recommendations

## 🎓 Educational Content

The documentation now serves as a comprehensive educational resource:

- **Theoretical background** for cavity-coupled molecular dynamics simulations
- **Mathematical formulations** for time-varying coupling
- **Physical interpretations** of simulation results
- **Best practices** for scientific computing
- **Performance optimization** strategies

This enhanced documentation provides a complete guide for researchers and developers working with cavity molecular dynamics simulations, from basic usage to advanced time-varying coupling experiments with comprehensive analysis capabilities. 