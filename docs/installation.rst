============
Installation
============

Requirements
============

- HOOMD-blue 4.0+
- Python 3.8+
- NumPy
- For GPU support: CUDA toolkit

Installation
============

.. code-block:: bash

   git clone https://github.com/muhammadhasyim/cav-hoomd.git
   cd cav-hoomd
   ./build_install.sh

**Installation Options**

With GPU support (default):

.. code-block:: bash

   ./build_install.sh

CPU only:

.. code-block:: bash

   ./build_install.sh --no-gpu

Uninstall
=========

To remove the installed plugins:

.. code-block:: bash

   ./uninstall.sh

Verify Installation
===================

.. code-block:: bash

   python examples/05_advanced_run.py --help

If you see the help message, the installation was successful!

**Test with Python:**

.. code-block:: python

   import hoomd
   from hoomd.cavitymd import CavityForce
   from hoomd.bussi_reservoir import BussiReservoir
   
   print("Cavity HOOMD installed successfully!")

Troubleshooting
===============

**Import errors:**
- Make sure HOOMD-blue is installed: ``python -c "import hoomd; print(hoomd.__version__)"``
- Check that the build completed without errors

**GPU issues:**
- Try CPU-only build: ``./build_install.sh --no-gpu``
- Verify CUDA installation: ``nvcc --version``

**Need help?**
- Check build logs for error messages
- Post an issue: https://github.com/muhammadhasyim/cav-hoomd/issues 