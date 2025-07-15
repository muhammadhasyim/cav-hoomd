============
Installation
============

System Requirements
===================

**Operating System:**
- Linux (Ubuntu 20.04+, RHEL 8+, or equivalent)
- macOS 10.15+ (experimental)

**Software Dependencies:**
- Python 3.8+
- HOOMD-blue 5.2.0 (see building instructions below)
- CUDA 11.0+ (optional, for GPU acceleration)
- CMake 3.15+
- C++ compiler with C++17 support (GCC 7+, Clang 5+)

**Python Dependencies:**
- NumPy 1.19+
- SciPy 1.5+
- matplotlib (for analysis)
- pandas (for data processing)
- gsd (for trajectory files)
- pybind11 (for C++ bindings)

Building HOOMD-blue from Source
===============================

**You need to build HOOMD-blue from source first** before installing the Cavity HOOMD plugin. Follow these commands:

**Install Prerequisites:**

.. code-block:: bash

   micromamba install cmake eigen git numpy pybind11 python

**Obtain the Source:**

.. code-block:: bash

   git clone --recursive git@github.com:glotzerlab/hoomd-blue.git --branch v5.2.0

**Change to the Repository Directory:**

.. code-block:: bash

   cd hoomd-blue

**Configure:**

.. code-block:: bash

   cmake -B build -S .

**Build the Package:**

.. code-block:: bash

   cd build
   make -j$(nproc)

**Run Tests:**

.. code-block:: bash

   python3 -m pytest hoomd

**Install the Package (optional):**

.. code-block:: bash

   make install

**Verify HOOMD-blue Version:**

.. code-block:: bash

   python3 -c "import hoomd; print('HOOMD version:', hoomd.version.version)"
   python3 -c "import hoomd; assert hoomd.version.version == '5.2.0', f'Expected HOOMD 5.2.0, got {hoomd.version.version}'"

.. note::
   After successfully building and installing HOOMD-blue from source, you can proceed with the Cavity HOOMD plugin installation below.

Quick Installation
==================

**Standard Installation (GPU + CPU):**

.. code-block:: bash

   git clone https://github.com/muhammadhasyim/cav-hoomd.git
   cd cav-hoomd
   ./build_install.sh

**CPU-only Installation:**

.. code-block:: bash

   git clone https://github.com/muhammadhasyim/cav-hoomd.git
   cd cav-hoomd
   ./build_install.sh --no-gpu

**Verify Installation:**

.. code-block:: bash

   python -c "import hoomd.cavitymd; print('Installation successful!')"
   python examples/05_advanced_run.py --help
   
   # Verify exact HOOMD-blue version compatibility
   python -c "import hoomd; assert hoomd.version.version == '5.2.0', f'ERROR: Expected HOOMD 5.2.0, got {hoomd.version.version}'"

If all commands work without errors, the installation is complete.

Detailed Installation Process
=============================

**Step 1: Clone Repository**

.. code-block:: bash

   git clone https://github.com/muhammadhasyim/cav-hoomd.git
   cd cav-hoomd

**Step 2: Choose Installation Method**

The simplified architecture supports automated installation:

.. code-block:: bash

   # Full installation with GPU support
   ./build_install.sh
   
   # CPU-only installation
   ./build_install.sh --no-gpu
   
   # Debug build for development
   ./build_install.sh --debug

**Step 3: Verify Components**

.. code-block:: bash

   # Test core plugin
   python -c "from hoomd.cavitymd.forces import CavityForce; print('CavityForce imported successfully')"
   
   # Test Bussi thermostat
   python -c "from hoomd.bussi_reservoir import BussiReservoir; print('BussiReservoir imported successfully')"
   
   # Test GPU support (if installed)
   python -c "import hoomd; print('GPU available:', hoomd.device.GPU.is_available())"

**Step 4: Run Test Simulation**

.. code-block:: bash

   python examples/05_advanced_run.py --coupling 1e-3 --runtime 100

Architecture Overview
=====================

**Simplified Plugin Structure**

The new architecture eliminates complex nested directories:

.. code-block:: text

   cav-hoomd/
   ├── src/
   │   ├── cavitymd/              # Main cavity MD plugin
   │   │   ├── __init__.py
   │   │   ├── forces.py          # CavityForce implementation
   │   │   ├── simulation.py      # CavityMDSimulation class
   │   │   ├── analysis.py        # Energy and correlation tracking
   │   │   ├── variants.py        # Time-varying coupling support
   │   │   └── ...
   │   ├── bussi_reservoir/       # Bussi thermostat plugin
   │   └── *.cc, *.h, *.cu       # C++/CUDA implementations
   ├── examples/
   │   ├── 05_advanced_run.py     # Main simulation script
   │   ├── 05_advanced_run.ipynb  # Jupyter notebook
   │   └── analysis scripts
   ├── docs/                      # Documentation
   └── build_install.sh           # Automated build script

**Key Improvements**

- **Direct imports**: ``from hoomd.cavitymd.forces import CavityForce``
- **Simplified building**: Single ``build_install.sh`` script
- **Better documentation**: Real docstrings preserved
- **Conda compatibility**: Works with Read the Docs and conda environments

Installation Options
====================

**Option 1: Automated Script (Recommended)**

.. code-block:: bash

   ./build_install.sh [--no-gpu] [--debug]

This script:
- Installs/updates HOOMD-blue
- Compiles C++/CUDA components
- Installs Python packages
- Verifies installation

**Option 2: Manual Installation**

.. code-block:: bash

   # Install HOOMD-blue
   conda install -c conda-forge hoomd
   
   # Build C++ components
   mkdir build && cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release
   make -j$(nproc)
   
   # Install Python packages
   pip install -e .

**Option 3: Development Installation**

.. code-block:: bash

   # Install in development mode
   ./build_install.sh --debug
   
   # This enables:
   # - Debug symbols
   # - Verbose output
   # - Development tools

GPU Support
===========

**CUDA Requirements:**
- CUDA 11.0+ (CUDA 12.0+ recommended)
- Compatible GPU (compute capability 6.0+)
- Sufficient GPU memory (>= 2GB recommended)

**GPU Installation:**

.. code-block:: bash

   # Check CUDA availability
   nvidia-smi
   nvcc --version
   
   # Install with GPU support
   ./build_install.sh

**GPU Testing:**

.. code-block:: bash

   # Test GPU functionality
   python examples/05_advanced_run.py --device GPU --coupling 1e-3 --runtime 100

**Performance Benefits:**

- **10-100x speedup** for large systems
- **Optimized kernels** for cavity force computation
- **Memory coalescing** for efficient GPU memory access
- **Automatic fallback** to CPU if GPU unavailable

High-Performance Computing
==========================

**SLURM Integration:**

The plugin includes native SLURM support for HPC environments:

.. code-block:: bash

   # Submit array job
   sbatch --array=1-10 submit.sh
   
   # Automatic replica detection
   python examples/05_advanced_run.py --coupling 1e-3 --runtime 1000

**Conda Environment Setup:**

.. code-block:: bash

   # Create isolated environment
   conda create -n cavity-hoomd python=3.9
   conda activate cavity-hoomd
   
   # Install in environment
   ./build_install.sh

**Module Loading (HPC Systems):**

.. code-block:: bash

   # Load required modules
   module load gcc/9.3.0
   module load cuda/11.4
   module load cmake/3.20.0
   
   # Build with specific compilers
   ./build_install.sh

Troubleshooting
===============

**Common Issues:**

1. **ImportError: No module named 'hoomd.cavitymd'**
   
   .. code-block:: bash
   
      # Reinstall from project root
      cd /path/to/cav-hoomd
      ./build_install.sh

2. **C++ compilation errors**
   
   .. code-block:: bash
   
      # Install build dependencies
      sudo apt-get install build-essential cmake
      
      # Try CPU-only build
      ./build_install.sh --no-gpu

3. **GPU not detected**
   
   .. code-block:: bash
   
      # Check CUDA installation
      nvidia-smi
      python -c "import hoomd; print(hoomd.device.GPU.is_available())"

4. **Permission errors**
   
   .. code-block:: bash
   
      # Use conda/virtual environment
      conda create -n cavity-hoomd python=3.9
      conda activate cavity-hoomd
      ./build_install.sh

**Getting Help:**

- **Documentation**: https://cav-hoomd.readthedocs.io/
- **Issues**: https://github.com/muhammadhasyim/cav-hoomd/issues
- **Discussions**: https://github.com/muhammadhasyim/cav-hoomd/discussions

**Development Installation:**

.. code-block:: bash

   # Clone with development tools
   git clone https://github.com/muhammadhasyim/cav-hoomd.git
   cd cav-hoomd
   
   # Install pre-commit hooks
   pip install pre-commit
   pre-commit install
   
   # Development build
   ./build_install.sh --debug

Uninstallation
==============

**Remove Installation:**

.. code-block:: bash

   # Run uninstall script
   ./uninstall.sh
   
   # Manual removal
   pip uninstall hoomd cavitymd bussi-reservoir

**Clean Build:**

.. code-block:: bash

   # Remove build artifacts
   rm -rf build/
   rm -rf __pycache__/
   rm -rf src/__pycache__/
   
   # Rebuild
   ./build_install.sh

Update Installation
===================

**Update to Latest Version:**

.. code-block:: bash

   # Pull latest changes
   git pull origin main
   
   # Rebuild and reinstall
   ./build_install.sh

**Version Information:**

.. code-block:: bash

   # Check installed version
   python -c "import hoomd.cavitymd; print('Cavity HOOMD installed')"
   
   # Check HOOMD-blue version (must be exactly 5.2.0)
   python -c "import hoomd; print('HOOMD version:', hoomd.version.version)"
   python -c "import hoomd; assert hoomd.version.version == '5.2.0', f'ERROR: Expected HOOMD 5.2.0, got {hoomd.version.version}'"

This simplified installation process ensures reliable builds across different systems while maintaining full functionality and performance optimization. 