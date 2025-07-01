================
Parameter Sweeps
================

**05_advanced_run.py** makes it easy to explore parameter spaces systematically. You can sweep over coupling strengths, temperatures, frequencies, and run multiple replicas with simple command-line options.

Basic Parameter Sweeps
======================

**Coupling Strength Sweep**

.. code-block:: bash

   # Sweep over multiple coupling strengths
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-4,1e-3,1e-2 --runtime 1000

This automatically creates separate directories for each coupling strength:

.. code-block:: text

   bussi_langevin_finiteq_coupling_1e-04/
   bussi_langevin_finiteq_coupling_1e-03/
   bussi_langevin_finiteq_coupling_1e-02/

**Temperature Sweep**

.. code-block:: bash

   # Sweep over multiple temperatures
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-3 --temperature 50,100,200,300 --runtime 1000

**Frequency Sweep**

.. code-block:: bash

   # Sweep over cavity frequencies
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-3 --frequency 1500,2000,2500,3000 --runtime 1000

Multi-Dimensional Sweeps
========================

**Combined Parameter Sweeps**

You can sweep over multiple parameters simultaneously:

.. code-block:: bash

   # Sweep coupling AND temperature
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-4,1e-3,1e-2 --temperature 100,200,300 --runtime 1000

This creates a directory for each combination:

.. code-block:: text

   bussi_langevin_finiteq_c1e-04_T100K_f2000cm/
   bussi_langevin_finiteq_c1e-04_T200K_f2000cm/
   bussi_langevin_finiteq_c1e-04_T300K_f2000cm/
   bussi_langevin_finiteq_c1e-03_T100K_f2000cm/
   ...

**Three-Parameter Sweep**

.. code-block:: bash

   # Sweep coupling, temperature, AND frequency
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-4,1e-3 --temperature 100,300 --frequency 2000,3000 --runtime 500

Multiple Replicas
=================

**Replica Ranges**

.. code-block:: bash

   # Run replicas 1-5 for each parameter combination
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-3,1e-2 --replicas "1-5" --runtime 1000

**Specific Replicas**

.. code-block:: bash

   # Run specific replica numbers
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-3 --replicas "1,3,5,7,9" --runtime 1000

**Large Replica Studies**

.. code-block:: bash

   # Run 20 replicas for good statistics
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-3 --replicas "1-20" --runtime 2000

HPC and SLURM Integration
=========================

**SLURM Array Jobs**

The script automatically detects SLURM environments and uses array task IDs:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --job-name=cavity_sweep
   #SBATCH --array=1-10
   #SBATCH --time=24:00:00
   #SBATCH --cpus-per-task=4

   # This automatically uses SLURM_ARRAY_TASK_ID as replica number
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-3,1e-4,1e-5 --runtime 5000

**Manual Replica Control**

If not using SLURM arrays, specify replicas manually:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --job-name=cavity_sweep
   #SBATCH --time=24:00:00

   # Run specific parameter combinations
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-3 --replicas "1-10" --runtime 5000

Common Sweep Patterns
=====================

**Coupling Strength Study**

.. code-block:: bash

   # Explore weak to strong coupling regime
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-5,1e-4,1e-3,1e-2,1e-1 --runtime 2000

**Temperature Dependence**

.. code-block:: bash

   # Study temperature effects
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-3 --temperature 50,75,100,150,200,300 --runtime 1500

**Resonance Study**

.. code-block:: bash

   # Scan around molecular resonance
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-3 --frequency 1800,1900,2000,2100,2200 --runtime 1500

**Thermostat Comparison**

.. code-block:: bash

   # Compare all thermostat combinations
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000
   python 05_advanced_run.py --experiment bussi_langevin_no_finiteq --coupling 1e-3 --runtime 1000
   python 05_advanced_run.py --experiment langevin_langevin --coupling 1e-3 --runtime 1000
   python 05_advanced_run.py --experiment bussi_bussi --coupling 1e-3 --runtime 1000

**Cavity vs Control**

.. code-block:: bash

   # Compare cavity vs no-cavity for multiple conditions
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-3 --temperature 100,200,300 --runtime 1000

   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --no-cavity --temperature 100,200,300 --runtime 1000

Advanced Sweep Options
======================

**Custom Output Periods**

.. code-block:: bash

   # Reduce output frequency for large sweeps
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-4,1e-3,1e-2 --runtime 5000 \
       --energy-output-period-ps 1.0 \
       --gsd-output-period-ps 100.0

**Analysis-Focused Sweeps**

.. code-block:: bash

   # Enable detailed tracking for parameter sweep
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-4,1e-3,1e-2 --runtime 2000 \
       --enable-energy-tracker --enable-fkt

**GPU Acceleration**

.. code-block:: bash

   # Use GPU for faster parameter sweeps
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-4,1e-3,1e-2 --runtime 3000 \
       --device GPU --gpu-id 0

Managing Large Sweeps
=====================

**Organizing Output**

For large parameter sweeps, output is automatically organized:

.. code-block:: text

   # Example structure for coupling + temperature sweep
   bussi_langevin_finiteq_c1e-04_T100K_f2000cm/
   ├── prod-1.gsd
   ├── prod-1-energy.txt
   └── prod-1.log
   
   bussi_langevin_finiteq_c1e-04_T200K_f2000cm/
   ├── prod-1.gsd
   ├── prod-1-energy.txt
   └── prod-1.log
   
   ...

**Storage Considerations**

.. code-block:: bash

   # Reduce storage for large sweeps
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-5,1e-4,1e-3,1e-2 --runtime 10000 \
       --gsd-output-period-ps 200.0  # Less frequent trajectory output

**Progress Monitoring**

.. code-block:: bash

   # Enable logging for sweep monitoring
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-4,1e-3,1e-2 --runtime 5000 \
       --log-to-file --log-to-console

Example Analysis Scripts
========================

**Collect Sweep Results**

.. code-block:: python

   import pandas as pd
   import glob
   import matplotlib.pyplot as plt

   # Collect energy data from all coupling strengths
   coupling_values = ['1e-04', '1e-03', '1e-02']
   all_data = {}

   for coupling in coupling_values:
       pattern = f"bussi_langevin_finiteq_coupling_{coupling}/prod-*-energy.txt"
       files = glob.glob(pattern)
       
       if files:
           data = pd.read_csv(files[0], delimiter='\t')
           all_data[coupling] = data

   # Plot energy drift vs coupling strength
   drifts = []
   for coupling, data in all_data.items():
       drift = (data['total_energy'].iloc[-1] - data['total_energy'].iloc[0]) / data['total_energy'].iloc[0]
       drifts.append((float(coupling.replace('e-0', 'e-')), abs(drift)))

   drifts.sort()
   couplings, drift_values = zip(*drifts)

   plt.figure(figsize=(10, 6))
   plt.loglog(couplings, drift_values, 'o-')
   plt.xlabel('Coupling Strength')
   plt.ylabel('|Energy Drift|')
   plt.title('Energy Conservation vs Coupling Strength')
   plt.show()

**Statistical Analysis**

.. code-block:: python

   # Analyze multiple replicas
   import numpy as np

   coupling = '1e-03'
   replica_data = []

   # Load all replicas for this coupling
   for replica in range(1, 6):  # Replicas 1-5
       pattern = f"bussi_langevin_finiteq_coupling_{coupling}/prod-{replica}-energy.txt"
       files = glob.glob(pattern)
       
       if files:
           data = pd.read_csv(files[0], delimiter='\t')
           final_energy = data['total_energy'].iloc[-1]
           replica_data.append(final_energy)

   # Calculate statistics
   mean_energy = np.mean(replica_data)
   std_energy = np.std(replica_data)
   
   print(f"Final energy: {mean_energy:.6f} ± {std_energy:.6f} Hartree")
   print(f"Relative std: {std_energy/abs(mean_energy)*100:.2f}%")

Best Practices
==============

**Sweep Planning**

1. **Start Small**: Test with short runtimes first
2. **Check One Parameter**: Verify one parameter works before combining
3. **Use Controls**: Always include no-cavity simulations
4. **Plan Storage**: Large sweeps generate lots of data

**Parameter Ranges**

- **Coupling**: Use logarithmic spacing (1e-5, 1e-4, 1e-3, ...)
- **Temperature**: Use linear or geometric spacing
- **Frequency**: Focus around molecular resonances

**Computational Efficiency**

.. code-block:: bash

   # Efficient sweep settings
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-4,1e-3,1e-2 --runtime 2000 \
       --gsd-output-period-ps 100.0 \
       --console-output-period-ps 10.0

**Error Handling**

- Monitor log files for failed simulations
- Use shorter runtimes for initial testing
- Check energy conservation for each parameter set

Next Steps
==========

* Learn about :doc:`analysis` to process your sweep results
* Check :doc:`troubleshooting` for handling sweep issues
* See :doc:`basic_usage` for more simulation options 