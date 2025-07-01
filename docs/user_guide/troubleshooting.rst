===============
Troubleshooting
===============

Common issues and solutions when using **05_advanced_run.py** for cavity MD simulations.

Installation Issues
===================

**"ModuleNotFoundError: No module named 'hoomd'"**

Install HOOMD-blue first:

.. code-block:: bash

   # Install HOOMD-blue
   conda install -c conda-forge hoomd

   # Then install cavity-hoomd
   pip install -e .

**"No module named 'hoomd.cavitymd'"**

The plugin isn't installed correctly:

.. code-block:: bash

   # Reinstall the plugin
   cd /path/to/cavity-hoomd
   pip install -e . --force-reinstall

**"ImportError: cannot import name 'CavityForce'"**

Check that the C++ extension compiled correctly:

.. code-block:: bash

   # Check if the plugin is available
   python -c "import hoomd; print(hoomd.version.plugins)"

File and Directory Issues
=========================

**"FileNotFoundError: [Errno 2] No such file or directory: '../init-0.gsd'"**

The script looks for initial structure files:

.. code-block:: bash

   # Check if init file exists
   ls init-0.gsd

   # Or create a symlink if file is elsewhere
   ln -s /path/to/your/initial.gsd init-0.gsd

**"ERROR: Cavity simulation requested but no cavity particle type 'L' found"**

Your GSD file doesn't contain a cavity particle:

.. code-block:: bash

   # Run with --no-cavity for molecular-only simulation
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --no-cavity --coupling 1e-3 --runtime 1000

   # Or the script will automatically add a cavity particle if needed

**"PermissionError: [Errno 13] Permission denied"**

Check directory permissions:

.. code-block:: bash

   # Make sure you can write to the current directory
   ls -la .

   # Or run in a different directory
   mkdir my_simulations
   cd my_simulations
   python ../05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000

Simulation Runtime Errors
=========================

**"RuntimeError: CUDA out of memory"**

GPU memory is exhausted:

.. code-block:: bash

   # Force CPU usage
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --device CPU --coupling 1e-3 --runtime 1000

   # Or use different GPU
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --device GPU --gpu-id 1 --coupling 1e-3 --runtime 1000

**"RuntimeError: Particle out of box"**

Particles have moved outside the simulation box:

.. code-block:: bash

   # Reduce coupling strength
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-4 --runtime 1000

   # Use fixed timestep with smaller steps
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --fixed-timestep --timestep 0.5 --coupling 1e-3 --runtime 1000

**"Simulation crashes without error message"**

Enable detailed logging:

.. code-block:: bash

   # Enable logging to see what's happening
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --log-to-file --log-to-console

**"ERROR: Cannot use Langevin thermostat with tau=0"**

Thermostat time constants must be positive:

.. code-block:: bash

   # Fix with proper tau values
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --molecular-tau 5.0 --cavity-tau 1.0

Performance Issues
==================

**Simulation runs very slowly**

Try these optimizations:

.. code-block:: bash

   # Use GPU if available
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --device GPU --coupling 1e-3 --runtime 1000

   # Reduce output frequency
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --gsd-output-period-ps 100.0 --console-output-period-ps 10.0

   # Use fixed timestep (can be faster than adaptive)
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --fixed-timestep --timestep 2.0 \
       --coupling 1e-3 --runtime 1000

**"ns/day performance is very low"**

Check system resources:

.. code-block:: bash

   # Monitor CPU/GPU usage
   htop    # For CPU
   nvidia-smi  # For GPU

   # Try different device
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --device CPU --coupling 1e-3 --runtime 1000

**"Memory usage keeps growing"**

Disable or reduce analysis features:

.. code-block:: bash

   # Disable energy tracker for long runs
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 10000 \
       --gsd-output-period-ps 200.0

   # Limit energy output time
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 10000 \
       --enable-energy-tracker --max-energy-output-time 1000

Energy Conservation Issues
=========================

**"Energy drift is too large"**

Check energy conservation:

.. code-block:: bash

   # Enable energy tracking to monitor
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --enable-energy-tracker

   # Use adaptive timestep (default) instead of fixed
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000

   # Reduce coupling strength
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-4 --runtime 1000

**"Total energy is not conserved"**

This might be expected with thermostats:

.. code-block:: bash

   # For energy conservation, use NVE (no thermostats)
   # Note: This requires manual setup, not available in 05_advanced_run.py presets

   # Check if reservoir energies are tracked properly
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --enable-energy-tracker

**"Kinetic energy fluctuates wildly"**

Temperature control issues:

.. code-block:: bash

   # Increase thermostat time constants for gentler coupling
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --molecular-tau 10.0 --cavity-tau 5.0

SLURM and HPC Issues
====================

**"SLURM array job not working correctly"**

Check SLURM environment:

.. code-block:: bash

   # Check if SLURM variables are set
   echo $SLURM_ARRAY_TASK_ID
   echo $SLURM_JOB_ID

   # Manually specify replica if needed
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --replicas "1"

**"Job runs out of time"**

Optimize for time limits:

.. code-block:: bash

   # Reduce runtime for testing
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 500

   # Reduce output frequency
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 2000 \
       --gsd-output-period-ps 200.0 --energy-output-period-ps 1.0

**"Multiple jobs writing to same directory"**

Make sure output directories are unique:

.. code-block:: bash

   # For parameter sweeps, directories are auto-generated
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3,1e-4 --runtime 1000

   # For replicas, use SLURM array task ID
   #SBATCH --array=1-10

Analysis Issues
===============

**"Output files are empty or corrupted"**

Check simulation completion:

.. code-block:: bash

   # Check log files for errors
   grep -i error */prod-*.log

   # Verify simulation completed
   tail */prod-*.log

**"Energy file missing columns"**

Some features may not be enabled:

.. code-block:: bash

   # Enable comprehensive energy tracking
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --enable-energy-tracker

**"F(k,t) file not generated"**

F(k,t) tracking must be explicitly enabled:

.. code-block:: bash

   # Enable F(k,t) correlation tracking
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --enable-fkt

**"Cannot read GSD file"**

Check GSD file integrity:

.. code-block:: python

   import gsd.hoomd
   
   try:
       with gsd.hoomd.open('prod-1.gsd', 'r') as f:
           print(f"Frames: {len(f)}")
   except:
       print("GSD file is corrupted")

Parameter Issues
================

**"No cavity effects observed"**

Check parameter values:

.. code-block:: bash

   # Increase coupling strength
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-2 --runtime 1000

   # Check cavity frequency matches molecular resonances
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --frequency 1800 --runtime 1000

   # Increase simulation time
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 5000

**"Simulation becomes unstable with strong coupling"**

Reduce timestep or coupling:

.. code-block:: bash

   # Use fixed smaller timestep
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-2 --runtime 1000 \
       --fixed-timestep --timestep 0.5

   # Start with weaker coupling
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000

**"Results don't match literature"**

Check parameter units and values:

.. code-block:: bash

   # Verify frequency is in cm⁻¹
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --frequency 2000 --runtime 1000

   # Check temperature in Kelvin
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --temperature 300 --runtime 1000

Command Line Issues
===================

**"Unknown experiment type"**

Check available experiments:

.. code-block:: bash

   # List available experiments
   python 05_advanced_run.py --help

   # Use correct experiment name
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000

**"Invalid coupling format"**

Use proper number format:

.. code-block:: bash

   # Correct formats
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 0.001 --runtime 1000
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3,1e-4 --runtime 1000

**"Replicas format error"**

Use proper replica specification:

.. code-block:: bash

   # Correct formats
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --replicas "1-5" --runtime 1000
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --replicas "1,3,5" --runtime 1000

Getting More Help
=================

**Enable Debug Output**

For detailed debugging:

.. code-block:: bash

   # Maximum logging output
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --log-to-file --log-to-console --enable-energy-tracker

**Check System Requirements**

Verify your environment:

.. code-block:: python

   import hoomd
   print(f"HOOMD version: {hoomd.version.version}")
   print(f"HOOMD plugins: {hoomd.version.plugins}")
   
   # Check for GPU support
   try:
       device = hoomd.device.GPU()
       print("GPU support available")
   except:
       print("GPU support not available")

**Contact Support**

If issues persist:

1. **GitHub Issues**: https://github.com/yourusername/cavity-hoomd/issues
2. **Include**: Full error messages, command used, system info
3. **Attach**: Log files and minimal reproducing example

**Useful Information to Include**

- Operating system and version
- Python version
- HOOMD-blue version
- Complete error message
- Command that caused the error
- Contents of log files

Next Steps
==========

* Review :doc:`basic_usage` for correct command syntax
* Check :doc:`parameter_sweeps` for parameter exploration
* See :doc:`analysis` for understanding output files 