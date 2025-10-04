#!/usr/bin/env python3
"""
Analysis of Harmonic Energy: Linear vs Power Law Temperature Dependence

This script demonstrates why harmonic oscillator energies should be linear in 
temperature and investigates the physical implications of the fitted models.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
import os

# Physical constants
k_B_hartree_per_K = 3.1668105e-6  # Boltzmann constant in Hartree/K

def linear_model(T, E0, slope):
    """Linear temperature dependence: E = E0 + slope*T"""
    return E0 + slope * T

def power_model(T, E0, A, alpha):
    """Power law temperature dependence: E = E0 + A*T^alpha"""
    return E0 + A * T**alpha

def quantum_harmonic_model(T, E0, N_eff, omega_eff):
    """
    Quantum harmonic oscillator model (simplified single frequency)
    E = E0 + N_eff * hbar*omega * [1/2 + 1/(exp(hbar*omega/kT) - 1)]
    """
    # Convert omega from K to energy units (assuming omega in K)
    beta_omega = omega_eff / T
    # Avoid overflow for large beta_omega
    beta_omega = np.clip(beta_omega, 0, 50)
    
    try:
        thermal_factor = 1.0 / (np.exp(beta_omega) - 1)
    except:
        thermal_factor = np.where(beta_omega > 10, 0, 1/beta_omega)  # High/low T limits
    
    return E0 + N_eff * k_B_hartree_per_K * omega_eff * (0.5 + thermal_factor)

def analyze_harmonic_energy(empirical_file):
    """Comprehensive analysis of harmonic energy temperature dependence."""
    
    print(" HARMONIC ENERGY TEMPERATURE DEPENDENCE ANALYSIS")
    print("=" * 60)
    
    # Load data
    data = pd.read_csv(empirical_file, sep=r'\s+', comment='#')
    T = data['temperature'].values
    E_harmonic = data['harmonic_hartree'].values
    
    print(f" Dataset: {len(T)} points, T = {T.min():.1f} - {T.max():.1f} K")
    print(f"    Harmonic energy: {E_harmonic.min():.6f} - {E_harmonic.max():.6f} Ha")
    
    # Fit different models
    print(f"\n FITTING DIFFERENT MODELS:")
    
    # 1. Linear model
    popt_linear, pcov_linear = curve_fit(linear_model, T, E_harmonic)
    E0_lin, slope_lin = popt_linear
    E_pred_linear = linear_model(T, *popt_linear)
    rmse_linear = np.sqrt(np.mean((E_harmonic - E_pred_linear)**2))
    r2_linear = 1 - np.var(E_harmonic - E_pred_linear) / np.var(E_harmonic)
    
    print(f"    LINEAR: E = {E0_lin:.6f} + {slope_lin:.6f}·T")
    print(f"      RMSE: {rmse_linear:.6f} Ha, R² = {r2_linear:.4f}")
    
    # 2. Power law model
    popt_power, pcov_power = curve_fit(power_model, T, E_harmonic, 
                                      p0=[E_harmonic.min(), 0.001, 1.0])
    E0_pow, A_pow, alpha_pow = popt_power
    E_pred_power = power_model(T, *popt_power)
    rmse_power = np.sqrt(np.mean((E_harmonic - E_pred_power)**2))
    r2_power = 1 - np.var(E_harmonic - E_pred_power) / np.var(E_harmonic)
    
    print(f"    POWER: E = {E0_pow:.6f} + {A_pow:.6f}·T^{alpha_pow:.3f}")
    print(f"      RMSE: {rmse_power:.6f} Ha, R² = {r2_power:.4f}")
    
    # 3. Quantum harmonic (attempt - may not converge)
    try:
        # Initial guess: reasonable phonon frequency ~300K, effective N modes
        popt_quantum, _ = curve_fit(quantum_harmonic_model, T, E_harmonic,
                                   p0=[E_harmonic.min(), 1000, 300],
                                   bounds=([E_harmonic.min()-0.01, 100, 50],
                                          [E_harmonic.min()+0.01, 5000, 1000]))
        E0_q, N_eff, omega_eff = popt_quantum
        E_pred_quantum = quantum_harmonic_model(T, *popt_quantum)
        rmse_quantum = np.sqrt(np.mean((E_harmonic - E_pred_quantum)**2))
        r2_quantum = 1 - np.var(E_harmonic - E_pred_quantum) / np.var(E_harmonic)
        
        print(f"    QUANTUM: N_eff = {N_eff:.0f}, ω_eff = {omega_eff:.0f} K")
        print(f"      RMSE: {rmse_quantum:.6f} Ha, R² = {r2_quantum:.4f}")
        quantum_fit = True
    except:
        print(f"    QUANTUM: Fit failed (data may be too limited)")
        quantum_fit = False
    
    # Physical analysis
    print(f"\n PHYSICAL ANALYSIS:")
    
    # Expected slope for classical harmonic oscillators
    N_atoms = 500  # Approximate
    N_modes = 3 * N_atoms  # 3N total modes
    expected_slope_classical = N_modes * k_B_hartree_per_K
    
    print(f"   Classical expectation (3N·k_B): {expected_slope_classical:.6f} Ha/K")
    print(f"   Fitted linear slope: {slope_lin:.6f} Ha/K")
    print(f"   Ratio (fitted/expected): {slope_lin/expected_slope_classical:.3f}")
    
    if slope_lin/expected_slope_classical < 0.5:
        print(f"     Slope much smaller than expected - indicates coupling/anharmonicity")
    elif slope_lin/expected_slope_classical > 0.8:
        print(f"    Slope close to classical expectation")
    
    # Power law analysis
    print(f"\n   Power law α = {alpha_pow:.3f}:")
    if alpha_pow < 0.9:
        print(f"     α < 1: Sublinear - UNPHYSICAL for harmonic oscillators")
    elif alpha_pow > 1.1:
        print(f"     α > 1: Superlinear - possible anharmonic effects")
    else:
        print(f"    α ≈ 1: Nearly linear - consistent with harmonic behavior")
    
    # Extrapolation comparison
    T_extrap = np.array([5, 10, 20, 50, 77])
    print(f"\n EXTRAPOLATION TO LOW TEMPERATURES:")
    print(f"T(K) | Linear(Ha) | Power(Ha) | Diff(Ha) | Diff(K)")
    print(f"-----|------------|-----------|----------|--------")
    
    for t in T_extrap:
        E_lin = linear_model(t, *popt_linear)
        E_pow = power_model(t, *popt_power)
        diff_energy = abs(E_lin - E_pow)
        # Convert energy difference to temperature difference using linear slope
        diff_temp = diff_energy / slope_lin if slope_lin > 0 else 0
        print(f"{t:4.0f} | {E_lin:10.6f} | {E_pow:9.6f} | {diff_energy:8.6f} | {diff_temp:6.1f}")
    
    # Create visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Data and fits
    T_smooth = np.linspace(T.min()-20, T.max()+20, 300)
    
    ax1.scatter(T, E_harmonic, color='blue', alpha=0.7, s=50, label='Empirical Data')
    ax1.plot(T_smooth, linear_model(T_smooth, *popt_linear), 'r-', 
             linewidth=2, label=f'Linear: E₀ + {slope_lin:.4f}·T')
    ax1.plot(T_smooth, power_model(T_smooth, *popt_power), 'g--', 
             linewidth=2, label=f'Power: E₀ + A·T^{alpha_pow:.3f}')
    
    # Highlight extrapolation region
    ax1.axvspan(0, T.min(), alpha=0.2, color='red', label='Extrapolation Zone')
    ax1.axvspan(5, 77, alpha=0.1, color='orange', label='Simulation Range')
    
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('Harmonic Energy (Ha)')
    ax1.set_title('Harmonic Energy vs Temperature: Model Comparison')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, T.max()+20)
    
    # Plot 2: Residuals
    ax2.scatter(T, E_harmonic - E_pred_linear, color='red', alpha=0.7, 
               label=f'Linear (RMSE={rmse_linear:.6f})')
    ax2.scatter(T, E_harmonic - E_pred_power, color='green', alpha=0.7, 
               label=f'Power (RMSE={rmse_power:.6f})')
    ax2.axhline(0, color='black', linestyle='-', alpha=0.5)
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('Residuals (Ha)')
    ax2.set_title('Model Residuals')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Low-T extrapolation detail
    T_low = np.linspace(5, T.min(), 100)
    E_low_linear = linear_model(T_low, *popt_linear)
    E_low_power = power_model(T_low, *popt_power)
    
    ax3.plot(T_low, E_low_linear, 'r-', linewidth=2, label='Linear Extrapolation')
    ax3.plot(T_low, E_low_power, 'g--', linewidth=2, label='Power Extrapolation')
    ax3.fill_between(T_low, E_low_linear, E_low_power, alpha=0.3, 
                     label='Extrapolation Uncertainty')
    
    ax3.set_xlabel('Temperature (K)')
    ax3.set_ylabel('Harmonic Energy (Ha)')
    ax3.set_title('Low-Temperature Extrapolation Comparison')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Temperature conversion (systemic temperature from energy)
    E_test = np.linspace(E_harmonic.min(), E_harmonic.max(), 100)
    T_from_E_linear = (E_test - E0_lin) / slope_lin
    T_from_E_power = ((E_test - E0_pow) / A_pow)**(1/alpha_pow)
    
    ax4.plot(E_test, T_from_E_linear, 'r-', linewidth=2, label='Linear Model')
    ax4.plot(E_test, T_from_E_power, 'g--', linewidth=2, label='Power Model')
    ax4.scatter(E_harmonic, T, color='blue', alpha=0.7, s=30, label='Data')
    
    ax4.set_xlabel('Harmonic Energy (Ha)')
    ax4.set_ylabel('Temperature (K)')
    ax4.set_title('Energy → Temperature Conversion')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('harmonic_linear_vs_power_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Summary and recommendations
    print(f"\n SUMMARY AND RECOMMENDATIONS:")
    print(f"   Statistical best fit: {'Power law' if rmse_power < rmse_linear else 'Linear'}")
    print(f"   Physical expectation: Linear (harmonic oscillators)")
    print(f"   Extrapolation difference: Up to {max(abs(linear_model(T_extrap, *popt_linear) - power_model(T_extrap, *popt_power))):.6f} Ha")
    
    if alpha_pow < 0.95:
        print(f"     RECOMMENDATION: Use linear model for physical consistency")
        print(f"      The sublinear power law (α={alpha_pow:.3f}) is unphysical for harmonic modes")
        print(f"      Consider: 1) anharmonic corrections, 2) mode coupling, 3) data artifacts")
    else:
        print(f"    Both models reasonable, power law α≈1 suggests nearly harmonic behavior")
    
    return {
        'linear_params': popt_linear,
        'power_params': popt_power,
        'rmse_linear': rmse_linear,
        'rmse_power': rmse_power,
        'alpha': alpha_pow
    }

if __name__ == "__main__":
    empirical_file = "final_nodiss_cavitymd/potential_energy_components_vs_temperature.txt"
    results = analyze_harmonic_energy(empirical_file)
