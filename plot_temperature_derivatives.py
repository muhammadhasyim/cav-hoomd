#!/usr/bin/env python3
"""
Plot temperature derivatives (dE/dT) for harmonic energy models.

This script demonstrates how different correction models behave in their
temperature derivatives, especially in the low-temperature limit approaching 0K.

Key physical constraint: dE/dT|_{T→0} = 0.25 × N × k_B = 0.000396 Ha/K

Author: Scientific Software Engineer
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import argparse
import warnings

# Physical constants
k_B_hartree_per_K = 3.1668115634556076e-06  # Boltzmann constant in Hartree/K
N = 500  # Number of particles
target_slope = 0.25 * N * k_B_hartree_per_K  # Target low-T slope

def load_harmonic_data(filename):
    """Load harmonic energy vs temperature data."""
    df = pd.read_csv(filename, delim_whitespace=True, comment='#')
    T = df['temperature'].values
    E = df['harmonic_hartree'].values
    return T, E

# Model functions (same as in plot_harmonic_energy_corrected.py)
def classical_linear(T):
    """E = 0.25*N*k_B*T"""
    return 0.25 * N * k_B_hartree_per_K * T

def power_law_constrained(T, A, alpha):
    """E = A*T^α, constrained to E(0)=0"""
    return A * (T**alpha)

def classical_plus_correction(T, alpha_corr, A_corr):
    """E = 0.25*N*k_B*T + A_corr*T^α_corr"""
    return 0.25 * N * k_B_hartree_per_K * T + A_corr * (T**alpha_corr)

def classical_plus_sqrt(T, A_sqrt):
    """E = 0.25*N*k_B*T + A*√T"""
    return 0.25 * N * k_B_hartree_per_K * T + A_sqrt * np.sqrt(T)

def classical_plus_T32(T, A_T32):
    """E = 0.25*N*k_B*T + A*T^(3/2)"""
    return 0.25 * N * k_B_hartree_per_K * T + A_T32 * (T**(1.5))

def tanh_constrained(T, T_trans, alpha_high):
    """Constrained tanh transition"""
    x = T / T_trans
    E_linear = 0.25 * N * k_B_hartree_per_K * T
    E_power = 0.25 * N * k_B_hartree_per_K * (T**alpha_high)
    
    weight = 0.5 * (1 + np.tanh(2 * (x - 1)))
    return (1 - weight) * E_linear + weight * E_power

def exponential_constrained(T, T_char):
    """E = 0.25*N*k_B*T*(1 - exp(-T/T_char))"""
    return 0.25 * N * k_B_hartree_per_K * T * (1 - np.exp(-T / T_char))

# Derivative functions
def d_classical_linear_dT(T):
    """dE/dT = 0.25*N*k_B (constant)"""
    return np.full_like(T, 0.25 * N * k_B_hartree_per_K)

def d_power_law_constrained_dT(T, A, alpha):
    """dE/dT = A*α*T^(α-1)"""
    return A * alpha * (T**(alpha - 1))

def d_classical_plus_correction_dT(T, alpha_corr, A_corr):
    """dE/dT = 0.25*N*k_B + A_corr*α_corr*T^(α_corr-1)"""
    return 0.25 * N * k_B_hartree_per_K + A_corr * alpha_corr * (T**(alpha_corr - 1))

def d_classical_plus_sqrt_dT(T, A_sqrt):
    """dE/dT = 0.25*N*k_B + A/(2√T)"""
    # Avoid division by zero at T=0
    T_safe = np.where(T > 0, T, np.finfo(float).eps)
    return 0.25 * N * k_B_hartree_per_K + A_sqrt / (2 * np.sqrt(T_safe))

def d_classical_plus_T32_dT(T, A_T32):
    """dE/dT = 0.25*N*k_B + (3/2)*A*T^(1/2)"""
    return 0.25 * N * k_B_hartree_per_K + 1.5 * A_T32 * np.sqrt(T)

def d_tanh_constrained_dT(T, T_trans, alpha_high):
    """Numerical derivative of tanh model"""
    dT = 1e-6
    return (tanh_constrained(T + dT, T_trans, alpha_high) - 
            tanh_constrained(T - dT, T_trans, alpha_high)) / (2 * dT)

def d_exponential_constrained_dT(T, T_char):
    """dE/dT = 0.25*N*k_B*(1 - exp(-T/T_char)) + 0.25*N*k_B*(T/T_char)*exp(-T/T_char)"""
    exp_term = np.exp(-T / T_char)
    return 0.25 * N * k_B_hartree_per_K * (1 - exp_term + (T / T_char) * exp_term)

def fit_models(T_data, E_data):
    """Fit all models to the data and return fitted parameters."""
    fitted_params = {}
    
    # Classical Linear (no fitting needed)
    fitted_params['Classical Linear'] = {}
    
    # Power Law
    try:
        popt, _ = curve_fit(power_law_constrained, T_data, E_data, 
                           bounds=([1e-8, 0.1], [1e-3, 3.0]), p0=[3e-4, 1.0])
        fitted_params['Power Law'] = {'A': popt[0], 'alpha': popt[1]}
    except:
        fitted_params['Power Law'] = {'A': 3e-4, 'alpha': 1.0}
    
    # Classical + Correction
    try:
        popt, _ = curve_fit(classical_plus_correction, T_data, E_data,
                           bounds=([0.1, -1e-4], [3.0, 1e-4]), p0=[1.3, 1e-6])
        fitted_params['Classical + Correction'] = {'alpha_corr': popt[0], 'A_corr': popt[1]}
    except:
        fitted_params['Classical + Correction'] = {'alpha_corr': 1.3, 'A_corr': 1e-6}
    
    # Classical + √T
    try:
        popt, _ = curve_fit(classical_plus_sqrt, T_data, E_data, p0=[1e-5])
        fitted_params['Classical + √T'] = {'A_sqrt': popt[0]}
    except:
        fitted_params['Classical + √T'] = {'A_sqrt': 1e-3}
    
    # Classical + T^(3/2)
    try:
        popt, _ = curve_fit(classical_plus_T32, T_data, E_data, p0=[1e-7])
        fitted_params['Classical + T^(3/2)'] = {'A_T32': popt[0]}
    except:
        fitted_params['Classical + T^(3/2)'] = {'A_T32': 1e-7}
    
    # Tanh Transition
    try:
        popt, _ = curve_fit(tanh_constrained, T_data, E_data,
                           bounds=([50, 0.1], [300, 3.0]), p0=[120, 0.8])
        fitted_params['Tanh Transition'] = {'T_trans': popt[0], 'alpha_high': popt[1]}
    except:
        fitted_params['Tanh Transition'] = {'T_trans': 120, 'alpha_high': 0.8}
    
    # Exponential
    try:
        popt, _ = curve_fit(exponential_constrained, T_data, E_data, p0=[100])
        fitted_params['Exponential'] = {'T_char': popt[0]}
    except:
        fitted_params['Exponential'] = {'T_char': 100}
    
    return fitted_params

def plot_derivatives(T_data, E_data, fitted_params, output_file=None, show_plot=True):
    """Plot temperature derivatives for all models."""
    
    # Create temperature range from 0K to max data temperature
    T_max = T_data.max()
    T_plot = np.linspace(0.1, T_max, 1000)  # Start from 0.1K to avoid singularities
    T_low = np.linspace(0.1, 20, 200)  # Focus on very low temperatures
    
    # Set up the plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    colors = ['red', 'blue', 'green', 'purple', 'magenta', 'orange', 'brown']
    linestyles = ['-', '-', '-', '-', '-', '--', '--']
    
    models_data = [
        ('Classical Linear', d_classical_linear_dT, fitted_params['Classical Linear'], {}),
        ('Power Law', d_power_law_constrained_dT, fitted_params['Power Law'], 
         {'A': fitted_params['Power Law']['A'], 'alpha': fitted_params['Power Law']['alpha']}),
        ('Classical + Correction', d_classical_plus_correction_dT, fitted_params['Classical + Correction'],
         {'alpha_corr': fitted_params['Classical + Correction']['alpha_corr'], 
          'A_corr': fitted_params['Classical + Correction']['A_corr']}),
        ('Classical + √T', d_classical_plus_sqrt_dT, fitted_params['Classical + √T'],
         {'A_sqrt': fitted_params['Classical + √T']['A_sqrt']}),
        ('Classical + T^(3/2)', d_classical_plus_T32_dT, fitted_params['Classical + T^(3/2)'],
         {'A_T32': fitted_params['Classical + T^(3/2)']['A_T32']}),
        ('Tanh Transition', d_tanh_constrained_dT, fitted_params['Tanh Transition'],
         {'T_trans': fitted_params['Tanh Transition']['T_trans'], 
          'alpha_high': fitted_params['Tanh Transition']['alpha_high']}),
        ('Exponential', d_exponential_constrained_dT, fitted_params['Exponential'],
         {'T_char': fitted_params['Exponential']['T_char']})
    ]
    
    # Plot 1: Full temperature range
    for i, (name, func, params, kwargs) in enumerate(models_data):
        try:
            if kwargs:
                dE_dT = func(T_plot, **kwargs)
            else:
                dE_dT = func(T_plot)
            
            ax1.plot(T_plot, dE_dT, color=colors[i], linestyle=linestyles[i], 
                    linewidth=2, label=name, alpha=0.8)
        except Exception as e:
            print(f"Warning: Could not plot {name}: {e}")
    
    # Add target slope line
    ax1.axhline(y=target_slope, color='black', linestyle=':', linewidth=2, 
               label=f'Target: {target_slope:.6f} Ha/K', alpha=0.7)
    
    ax1.set_xlabel('Temperature (K)', fontsize=12)
    ax1.set_ylabel('dE/dT (Ha/K)', fontsize=12)
    ax1.set_title('Temperature Derivatives: Full Range', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, T_max)
    
    # Add data range indicator
    ax1.axvspan(T_data.min(), T_data.max(), alpha=0.1, color='gray', label='Data Range')
    
    # Plot 2: Low temperature focus (0-20K)
    for i, (name, func, params, kwargs) in enumerate(models_data):
        try:
            if kwargs:
                dE_dT_low = func(T_low, **kwargs)
            else:
                dE_dT_low = func(T_low)
            
            # Handle infinite values for √T model
            if name == 'Classical + √T':
                # Clip extremely high values for visualization
                dE_dT_low = np.clip(dE_dT_low, 0, 0.01)
            
            ax2.plot(T_low, dE_dT_low, color=colors[i], linestyle=linestyles[i], 
                    linewidth=2, label=name, alpha=0.8)
        except Exception as e:
            print(f"Warning: Could not plot {name} at low T: {e}")
    
    # Add target slope line
    ax2.axhline(y=target_slope, color='black', linestyle=':', linewidth=2, 
               label=f'Target: {target_slope:.6f} Ha/K', alpha=0.7)
    
    ax2.set_xlabel('Temperature (K)', fontsize=12)
    ax2.set_ylabel('dE/dT (Ha/K)', fontsize=12)
    ax2.set_title('Temperature Derivatives: Low Temperature Focus (0-20K)', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 20)
    
    # Add extrapolation region indicator
    ax2.axvspan(0, T_data.min(), alpha=0.15, color='red', label='Deep Extrapolation')
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"💾 Plot saved to: {output_file}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()

def analyze_low_temperature_behavior(fitted_params):
    """Analyze behavior at very low temperatures."""
    print("\n🔍 LOW-TEMPERATURE DERIVATIVE ANALYSIS")
    print("=" * 60)
    print(f"🎯 Target: dE/dT|_{{T→0}} = {target_slope:.6f} Ha/K")
    print()
    
    T_test = np.array([0.1, 1.0, 5.0, 10.0])  # Test temperatures
    
    models_data = [
        ('Classical Linear', d_classical_linear_dT, fitted_params['Classical Linear'], {}),
        ('Power Law', d_power_law_constrained_dT, fitted_params['Power Law'], 
         {'A': fitted_params['Power Law']['A'], 'alpha': fitted_params['Power Law']['alpha']}),
        ('Classical + Correction', d_classical_plus_correction_dT, fitted_params['Classical + Correction'],
         {'alpha_corr': fitted_params['Classical + Correction']['alpha_corr'], 
          'A_corr': fitted_params['Classical + Correction']['A_corr']}),
        ('Classical + √T', d_classical_plus_sqrt_dT, fitted_params['Classical + √T'],
         {'A_sqrt': fitted_params['Classical + √T']['A_sqrt']}),
        ('Classical + T^(3/2)', d_classical_plus_T32_dT, fitted_params['Classical + T^(3/2)'],
         {'A_T32': fitted_params['Classical + T^(3/2)']['A_T32']}),
        ('Tanh Transition', d_tanh_constrained_dT, fitted_params['Tanh Transition'],
         {'T_trans': fitted_params['Tanh Transition']['T_trans'], 
          'alpha_high': fitted_params['Tanh Transition']['alpha_high']}),
        ('Exponential', d_exponential_constrained_dT, fitted_params['Exponential'],
         {'T_char': fitted_params['Exponential']['T_char']})
    ]
    
    for name, func, params, kwargs in models_data:
        print(f"\n📊 {name}:")
        try:
            if kwargs:
                derivatives = func(T_test, **kwargs)
            else:
                derivatives = func(T_test)
            
            for i, T in enumerate(T_test):
                ratio = (derivatives[i] / target_slope) * 100
                status = "✅" if 80 <= ratio <= 120 else "⚠️" if 50 <= ratio <= 150 else "❌"
                print(f"   T={T:4.1f}K: dE/dT = {derivatives[i]:.6f} Ha/K ({ratio:5.1f}%) {status}")
            
            # Check limit behavior
            limit_behavior = "→ 0.25×N×k_B" if name in ['Classical Linear', 'Classical + T^(3/2)', 'Tanh Transition'] else \
                           "→ ∞" if name == 'Classical + √T' else \
                           "→ A×α×T^(α-1)" if name == 'Power Law' else \
                           "→ varies"
            print(f"   Limit as T→0: {limit_behavior}")
            
        except Exception as e:
            print(f"   Error: {e}")

def main():
    parser = argparse.ArgumentParser(description='Plot temperature derivatives for harmonic energy models')
    parser.add_argument('--data-file', default='final_nodiss_cavitymd/potential_energy_components_vs_temperature.txt',
                       help='Input data file')
    parser.add_argument('--output', default='temperature_derivatives.png',
                       help='Output plot filename')
    parser.add_argument('--no-show', action='store_true',
                       help='Do not display the plot')
    
    args = parser.parse_args()
    
    print("🌡️ TEMPERATURE DERIVATIVE ANALYSIS")
    print("=" * 50)
    
    # Load data
    T_data, E_data = load_harmonic_data(args.data_file)
    print(f"📊 Loaded {len(T_data)} data points")
    print(f"   Temperature range: {T_data.min():.1f} - {T_data.max():.1f} K")
    print(f"🎯 Target slope: {target_slope:.6f} Ha/K")
    
    # Fit models
    print("\n🔧 Fitting models...")
    fitted_params = fit_models(T_data, E_data)
    
    # Plot derivatives
    plot_derivatives(T_data, E_data, fitted_params, 
                    output_file=args.output, show_plot=not args.no_show)
    
    # Analyze low-temperature behavior
    analyze_low_temperature_behavior(fitted_params)
    
    print("\n✅ Derivative analysis completed!")

if __name__ == "__main__":
    main()
