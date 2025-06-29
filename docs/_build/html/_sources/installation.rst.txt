============
Installation
============

System Requirements
===================

Cavity HOOMD requires:

* Python 3.8 or later
* HOOMD-blue 4.0 or later
* NumPy 1.19 or later
* SciPy 1.6 or later

Optional dependencies:

* CUDA 11.0+ (for GPU acceleration)
* MPI (for parallel simulations)

Installation Methods
====================

From Source (Recommended)
--------------------------

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/muhammadhasyim/cav-hoomd.git
   cd cav-hoomd

   # Install in development mode
   pip install -e .

   # Build and install the plugin
   ./build_install.sh

Verify Installation
===================

.. code-block:: python

   import hoomd
   from hoomd.cavitymd import CavityForce, PhysicalConstants

   print("Cavity HOOMD installed successfully!")
   print(f"Physical constants loaded: {PhysicalConstants.KB_HARTREE_PER_K}")

Troubleshooting
===============

Common Issues
-------------

**ImportError: No module named 'hoomd.cavitymd'**
   - Ensure HOOMD-blue is properly installed
   - Check that the plugin was built and installed correctly
   - Verify your Python environment

**CUDA errors**
   - Check CUDA installation and compatibility
   - Verify GPU drivers are up to date
   - Try CPU-only installation first

**Build errors**
   - Ensure all dependencies are installed
   - Check compiler compatibility
   - See the build logs for specific error messages

Getting Help
============

If you encounter issues:

1. Check the troubleshooting section above
2. Search existing issues on GitHub
3. Create a new issue with:
   - Your system information
   - Complete error messages
   - Minimal example that reproduces the problem 