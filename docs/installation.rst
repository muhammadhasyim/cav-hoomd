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

   git clone https://github.com/yourusername/cavity-hoomd.git
   cd cavity-hoomd
   cmake -B build -S .
   cmake --build build
   cmake --install build

**Installation Options**

With GPU support (default):

.. code-block:: bash

   cmake -B build -S . -DENABLE_GPU=ON
   cmake --build build
   cmake --install build

CPU only:

.. code-block:: bash

   cmake -B build -S . -DENABLE_GPU=OFF
   cmake --build build
   cmake --install build

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
- Try CPU-only build: ``cmake -B build -S . -DENABLE_GPU=OFF``
- Verify CUDA installation: ``nvcc --version``

**Need help?**
- Check build logs for error messages
- Post an issue: https://github.com/yourusername/cavity-hoomd/issues 