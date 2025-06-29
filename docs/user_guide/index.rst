==========
User Guide
==========

This comprehensive user guide will walk you through using Cavity HOOMD for cavity-coupled molecular dynamics simulations.

.. toctree::
   :maxdepth: 2

   basic_usage
   advanced_features
   parameter_sweeps
   analysis
   performance

Overview
========

Cavity HOOMD provides a complete framework for simulating molecular systems coupled to optical cavity modes. This guide covers:

**Getting Started**
   Learn the basic concepts and run your first simulation

**Advanced Features**
   Explore sophisticated capabilities like adaptive timesteps and energy tracking

**Parameter Sweeps**
   Efficiently explore parameter spaces with automated sweeps

**Analysis Tools**
   Understand and analyze simulation results

**Performance Optimization**
   Get the best performance from your simulations

Quick Navigation
================

.. grid:: 2
    :gutter: 3

    .. grid-item-card:: 🎯 **Basic Usage**
        :link: basic_usage
        :link-type: doc

        Start here to learn fundamental concepts and run your first
        cavity-coupled simulation.

    .. grid-item-card:: ⚡ **Advanced Features**
        :link: advanced_features
        :link-type: doc

        Explore adaptive timesteps, energy tracking, and other
        sophisticated simulation capabilities.

    .. grid-item-card:: 🔄 **Parameter Sweeps**
        :link: parameter_sweeps
        :link-type: doc

        Learn how to efficiently explore parameter spaces with
        automated sweeps and parallel execution.

    .. grid-item-card:: 📊 **Analysis Tools**
        :link: analysis
        :link-type: doc

        Master the analysis capabilities including energy tracking,
        F(k,t) correlations, and cavity mode monitoring.

Learning Path
=============

If you're new to Cavity HOOMD, we recommend following this learning path:

1. **Start with Basic Usage** (:doc:`basic_usage`)
   - Understand cavity-molecule coupling
   - Set up your first simulation
   - Interpret basic results

2. **Explore Advanced Features** (:doc:`advanced_features`)
   - Enable adaptive timestepping
   - Configure comprehensive energy tracking
   - Use different thermostat combinations

3. **Learn Parameter Sweeps** (:doc:`parameter_sweeps`)
   - Automate parameter exploration
   - Handle multiple replicas
   - Optimize computational resources

4. **Master Analysis Tools** (:doc:`analysis`)
   - Monitor energy conservation
   - Calculate correlation functions
   - Analyze cavity mode dynamics

5. **Optimize Performance** (:doc:`performance`)
   - Choose appropriate hardware
   - Configure for maximum efficiency
   - Troubleshoot common issues

Key Concepts
============

Before diving into the detailed guides, here are the key concepts you should understand:

**Cavity Modes**
   Optical cavity modes are represented as harmonic oscillators that can exchange energy with molecular degrees of freedom.

**Coupling Strength**
   The parameter `g` that determines how strongly the molecules interact with the cavity field.

**Finite vs Infinite Q**
   Whether the cavity mode can have finite momentum (finite-q) or is constrained to zero momentum (q=0).

**Thermostats**
   Independent temperature control for molecular and cavity degrees of freedom using various algorithms.

**Energy Conservation**
   The total energy (molecular + cavity + reservoir) should be conserved in properly configured simulations.

**Autocorrelation Functions**
   Time-dependent correlation functions like F(k,t) that reveal dynamic properties of the system.

Common Workflows
================

Here are some typical workflows you might follow:

**Basic Cavity Coupling Study**

1. Set up system with cavity coupling enabled
2. Run simulation with energy tracking
3. Analyze energy conservation and cavity mode dynamics
4. Compare with no-cavity control simulation

**Parameter Exploration**

1. Define parameter ranges (coupling, frequency, temperature)
2. Set up automated parameter sweep
3. Run multiple replicas in parallel
4. Analyze trends and identify optimal parameters

**Correlation Function Analysis**

1. Enable F(k,t) tracking during simulation
2. Run for sufficient time to capture relaxation
3. Post-process correlation data
4. Fit exponential decays and extract timescales

**Performance Optimization**

1. Test different hardware configurations
2. Optimize timestep and output frequencies
3. Profile memory usage and I/O patterns
4. Scale to larger systems or longer times

Getting Help
============

If you encounter issues or have questions:

1. **Check the FAQ** - Common issues and solutions
2. **Review Examples** - Working code for typical scenarios  
3. **API Reference** - Detailed parameter descriptions
4. **GitHub Issues** - Report bugs or request features
5. **Discussions** - Ask questions and share experiences

Next Steps
==========

Ready to get started? Head to :doc:`basic_usage` to run your first cavity-coupled simulation! 