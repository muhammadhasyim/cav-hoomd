#!/usr/bin/env python3
"""
Analyze Rosenfeld-Tarazona Extrapolation to Low Temperatures

This script analyzes how well the empirical energy-temperature relationship
extrapolates to low temperatures (5K-65K) for cavity MD simulations.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path

# Add the source directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

try:
    from cavitymd.analysis import EmpiricalTemperatureData
    EMPIRICAL_DATA_AVAILABLE = True
except ImportError:
    print("Warning: Could not import EmpiricalTemperatureData")
    EMPIRICAL_DATA_AVAILABLE = False

def analyze_extrapolation(empirical_data_file, energy_component='harmonic', output_file=None):
    """
    Analyze Rosenfeld-Tarazona extrapolation behavior.
    
    Parameters
    ----------
    empirical_data_file : str
        Path to empirical energy-temperature data file
    energy_component : str
        Energy component to analyze ('harmonic', 'lj_coulombic', 'total_PE')
    output_file : str, optional
        Output file for saving the plot
    """
    
    if not EMPIRICAL_DATA_AVAILABLE:
        print(" EmpiricalTemperatureData not available. Cannot perform analysis.")
        return
    
    print(f" Analyzing {energy_component} extrapolation...")
    
    # Load empirical data
    try:
        empirical_data = EmpiricalTemperatureData(
            data_file_path=empirical_data_file,
            energy_component=energy_component,
            use_direct_harmonic=False  # Use fitted Rosenfeld function
        )
    except Exception as e:
        print(f" Failed to load empirical data: {e}")
        return
    
    if not empirical_data.is_fitted:
        print(" Rosenfeld fitting failed. Cannot perform extrapolation analysis.")
        return
    
    # Extract fit parameters
    E0 = empirical_data.fit_params['E0']
    A = empirical_data.fit_params['A']
    alpha = empirical_data.fit_params['alpha']
    
    print(f"\n Rosenfeld-Tarazona Fit Parameters ({energy_component}):")
    print(f"   E₀ = {E0:.6f} Ha (0K energy)")
    print(f"   A  = {A:.6e} Ha/K^α")
    print(f"   α  = {alpha:.3f}")
    
    # Create temperature ranges
    T_empirical = empirical_data.temperatures
    T_min_empirical = T_empirical.min()
    T_max_empirical = T_empirical.max()
    
    # Extended temperature range for extrapolation
    T_extended = np.linspace(0, T_max_empirical + 50, 1000)
    T_extrapolation = T_extended[T_extended < T_min_empirical]
    T_simulation_range = np.linspace(5, 77, 100)  # Your actual simulation range
    
    # Calculate energies using Rosenfeld function
    def rosenfeld_function(T):
        return E0 + A * T**alpha
    
    E_extended = rosenfeld_function(T_extended)
    E_empirical_fit = rosenfeld_function(T_empirical)
    
    # Calculate temperatures from energies (inverse function)
    E_test_range = np.linspace(E0, empirical_data.energies.max(), 1000)
    T_from_energy = np.zeros_like(E_test_range)
    
    for i, E in enumerate(E_test_range):
        T_from_energy[i] = empirical_data.calculate_systemic_temperature(E)
    
    # Create the analysis plot
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Energy vs Temperature with extrapolation
    ax1.scatter(T_empirical, empirical_data.energies, 
               color='blue', alpha=0.7, s=50, label='Empirical Data')
    ax1.plot(T_extended, E_extended, 'r-', linewidth=2, 
             label=f'RT Fit: E₀+A·T^{alpha:.3f}')
    ax1.axvline(T_min_empirical, color='orange', linestyle='--', alpha=0.8,
               label=f'Data Range: {T_min_empirical:.0f}K')
    ax1.axvspan(5, T_min_empirical, alpha=0.2, color='red', 
               label='Extrapolation Zone')
    ax1.axvspan(5, 77, alpha=0.1, color='green', 
               label='Simulation Range')
    
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel(f'{energy_component.title()} Energy (Ha)')
    ax1.set_title('Energy vs Temperature: Empirical Data & Extrapolation')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, T_max_empirical + 20)
    
    # Plot 2: Temperature vs Energy (inverse relationship)
    ax2.plot(E_test_range, T_from_energy, 'g-', linewidth=2,
             label='T from E (Inverse RT)')
    ax2.scatter(empirical_data.energies, T_empirical, 
               color='blue', alpha=0.7, s=50, label='Empirical Data')
    ax2.axhspan(5, T_min_empirical, alpha=0.2, color='red', 
               label='Extrapolation Zone')
    ax2.axhspan(5, 77, alpha=0.1, color='green', 
               label='Simulation Range')
    
    ax2.set_xlabel(f'{energy_component.title()} Energy (Ha)')
    ax2.set_ylabel('Temperature (K)')
    ax2.set_title('Temperature vs Energy: Inverse Extrapolation')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Extrapolation validation at 77K
    if 77 < T_min_empirical:  # 77K is in extrapolation range
        T_77K = 77.0
        E_77K_predicted = rosenfeld_function(T_77K)
        T_77K_from_energy = empirical_data.calculate_systemic_temperature(E_77K_predicted)
        
        ax3.scatter([T_77K], [E_77K_predicted], color='red', s=100, 
                   label=f'77K Prediction: {E_77K_predicted:.6f} Ha')
        ax3.plot(T_extended, E_extended, 'r-', alpha=0.5)
        ax3.scatter(T_empirical, empirical_data.energies, 
                   color='blue', alpha=0.7, s=30)
        ax3.axvline(T_min_empirical, color='orange', linestyle='--')
        ax3.set_xlim(65, 85)
        ax3.set_ylim(E_77K_predicted - 0.01, empirical_data.energies.min() + 0.01)
        
        print(f"\n 77K Extrapolation Validation:")
        print(f"   Predicted Energy at 77K: {E_77K_predicted:.6f} Ha")
        print(f"   Round-trip Temperature: {T_77K_from_energy:.1f} K")
        print(f"   Error: {abs(T_77K_from_energy - T_77K):.2f} K")
    else:
        ax3.text(0.5, 0.5, '77K within empirical range', 
                transform=ax3.transAxes, ha='center', va='center', fontsize=14)
    
    ax3.set_xlabel('Temperature (K)')
    ax3.set_ylabel(f'{energy_component.title()} Energy (Ha)')
    ax3.set_title('77K Extrapolation Detail')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Residuals and fit quality
    residuals = empirical_data.energies - E_empirical_fit
    ax4.scatter(T_empirical, residuals, color='purple', alpha=0.7, s=50)
    ax4.axhline(0, color='black', linestyle='-', alpha=0.5)
    ax4.fill_between([T_min_empirical-10, T_min_empirical], 
                     [residuals.min()*1.5, residuals.min()*1.5],
                     [residuals.max()*1.5, residuals.max()*1.5],
                     alpha=0.2, color='red', label='Extrapolation Zone')
    
    rmse = np.sqrt(np.mean(residuals**2))
    ax4.set_xlabel('Temperature (K)')
    ax4.set_ylabel('Residuals (Ha)')
    ax4.set_title(f'Fit Quality: RMSE = {rmse:.6f} Ha')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    print(f"\n Fit Quality Statistics:")
    print(f"   RMSE: {rmse:.6f} Ha")
    print(f"   Max Residual: {np.abs(residuals).max():.6f} Ha")
    print(f"   R²: {1 - np.var(residuals)/np.var(empirical_data.energies):.4f}")
    
    # Physical validation checks
    print(f"\n Physical Validation:")
    
    # Check if E0 is reasonable (should be the lowest energy)
    E_min_empirical = empirical_data.energies.min()
    print(f"   E₀ vs Min Data: {E0:.6f} vs {E_min_empirical:.6f} Ha")
    if E0 > E_min_empirical:
        print("     Warning: E₀ > minimum data energy")
    else:
        print("    E₀ < minimum data energy (good)")
    
    # Check alpha value
    if 0.3 < alpha < 2.0:
        print(f"    α = {alpha:.3f} in reasonable range (0.3-2.0)")
    else:
        print(f"     α = {alpha:.3f} outside typical range")
    
    # Extrapolation range assessment
    extrap_range = T_min_empirical - 5  # 5K is your minimum simulation temperature
    print(f"   Extrapolation range: {extrap_range:.0f} K")
    if extrap_range > 50:
        print("     Large extrapolation range (>50K)")
    else:
        print("    Moderate extrapolation range")
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\n Plot saved to: {output_file}")
    
    plt.show()
    
    return {
        'fit_params': empirical_data.fit_params,
        'rmse': rmse,
        'extrapolation_range': extrap_range,
        'empirical_data': empirical_data
    }

def main():
    # Default empirical data file
    empirical_file = "final_nodiss_cavitymd/potential_energy_components_vs_temperature.txt"
    
    print(" ROSENFELD-TARAZONA EXTRAPOLATION ANALYSIS")
    print("=" * 60)
    
    # Analyze different energy components
    components = ['harmonic', 'lj_coulombic']
    
    for component in components:
        print(f"\n{'='*20} {component.upper()} {'='*20}")
        try:
            results = analyze_extrapolation(
                empirical_data_file=empirical_file,
                energy_component=component,
                output_file=f"extrapolation_analysis_{component}.png"
            )
        except Exception as e:
            print(f" Analysis failed for {component}: {e}")
        
        print("\n" + "-"*60)

if __name__ == "__main__":
    main()
