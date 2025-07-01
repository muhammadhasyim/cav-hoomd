======
Theory
======

Cavity HOOMD simulates molecules interacting with optical cavity modes. This enables study of light-matter interactions where molecular vibrations couple strongly to photons.

What It Models
==============

**Cavity-Molecule Interaction**

The simulation implements the interaction Hamiltonian:

.. math::

   H = \frac{1}{2} K q^2 + g \vec{q} \cdot \vec{d} + \frac{g^2}{2K} d^2

Where:
- :math:`q` is the cavity mode position (photon displacement)
- :math:`d` is the molecular dipole moment  
- :math:`g` is the coupling strength
- :math:`K = m \omega_c^2` is the cavity spring constant

**Physical Meaning**

- **First term**: Cavity harmonic oscillator energy
- **Second term**: Direct coupling between cavity and molecular dipole
- **Third term**: Dipole self-energy (A²/2 term from QED)

How It Works
============

**Classical Treatment**

The cavity mode is treated as a classical harmonic oscillator with:
- Position :math:`q` (photon displacement)
- Momentum :math:`p` (photon momentum)
- Natural frequency :math:`\omega_c`

**Forces**

Molecular particles feel additional forces:

.. math::

   F_i = -g \, q_i \, q - \frac{g^2}{K} q_i \vec{d}

The cavity particle experiences:

.. math::

   F_{\text{cavity}} = -K q - g \vec{d}

Key Physics
===========

**Strong Coupling Regime**

When :math:`g > \sqrt{\gamma \kappa}` (coupling stronger than loss rates), the system enters the strong coupling regime where:

- Molecular vibrations and cavity photons hybridize  
- New energy levels appear (polaritons)
- Collective effects emerge from many molecules

**Energy Conservation**

The total energy is conserved:

.. math::

   E_{\text{total}} = E_{\text{molecular}} + E_{\text{cavity}} + E_{\text{coupling}}

**Thermostats**

Different thermostat combinations control temperature:
- **Bussi**: Deterministic velocity rescaling
- **Langevin**: Stochastic friction and random forces

Applications
============

**What You Can Study**

- How cavity coupling affects molecular dynamics
- Energy transfer between molecules and light
- Collective vibrational effects
- Modified chemical reaction rates
- Polariton formation and dynamics

**Typical Parameter Ranges**

- Coupling strength: :math:`10^{-5}` to :math:`10^{-2}` (atomic units)
- Cavity frequency: 1000-3000 cm⁻¹ (molecular vibrations)
- Temperature: 50-300 K

**Observables**

The simulation tracks:
- Individual energy components (harmonic, coupling, dipole)
- Cavity mode position and momentum
- Molecular trajectory and dipole moment
- Total energy conservation

This framework enables computational study of cavity quantum electrodynamics effects in realistic molecular systems. 