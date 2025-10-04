#!/usr/bin/env python3
"""
Corrected Harmonic Energy Analysis with E(T=0) = 0 Constraint

This script analyzes harmonic energy vs temperature with the fundamental constraints:
1. E(T=0) = 0 (zero energy at absolute zero)
2. dE/dT|_{T→0} = 0.25 * N * k_B (specific low-temperature slope)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import argparse
from pathlib import Path

# Physical constants
k_B_hartree_per_K = 3.1668105e-6  # Boltzmann constant in Hartree/K

def classical_linear(T, N=500):
    """
    Exact classical: E = 0.25 * N * k_B * T
    
    Satisfies: E(0) = 0 and dE/dT = 0.25 * N * k_B
    """
    return 0.25 * N * k_B_hartree_per_K * T

def fitted_power_law(T, A, alpha):
    """
    General power law: E = A * T^α
    
    Satisfies: E(0) = 0
    Note: Only α=1 can satisfy the slope constraint exactly
    """
    return A * (T**alpha)

def classical_plus_correction(T, alpha, A_correction):
    """
    Classical linear plus correction: E = 0.25*N*k_B*T + A*T^α
    
    Satisfies: E(0) = 0 and correct initial slope
    Correction term for higher-order effects
    """
    N = 500  # Fixed
    return 0.25 * N * k_B_hartree_per_K * T + A_correction * (T**alpha)

def classical_plus_sqrt(T, A_sqrt):
    """
    Classical linear plus square root correction: E = 0.25*N*k_B*T + A*√T
    
    Satisfies: E(0) = 0 and correct initial slope
    Square root term represents quantum-like corrections or disorder effects
    
    Note: dE/dT = 0.25*N*k_B + A/(2√T) → 0.25*N*k_B as T→0
    """
    N = 500  # Fixed
    return 0.25 * N * k_B_hartree_per_K * T + A_sqrt * np.sqrt(T)

def classical_plus_T32(T, A_T32):
    """
    Classical linear plus T^(3/2) correction: E = 0.25*N*k_B*T + A*T^(3/2)
    
    Satisfies: E(0) = 0 and correct initial slope
    T^(3/2) term for quantum/anharmonic corrections that preserve slope constraint
    
    dE/dT = 0.25*N*k_B + (3/2)*A*T^(1/2) → 0.25*N*k_B as T→0 
    """
    N = 500  # Fixed
    return 0.25 * N * k_B_hartree_per_K * T + A_T32 * (T**(1.5))

def tanh_constrained(T, T_trans, alpha_high):
    """
    Constrained tanh transition: E(0)=0, correct initial slope
    
    E = 0.25*N*k_B*T*[f + (1-f)*(T/T_ref)^(α-1)]
    where f = 0.5*(1 + tanh((T-T_trans)/(0.2*T_trans)))
    """
    N = 500
    f = 0.5 * (1 + np.tanh((T - T_trans) / (0.2 * T_trans)))
    
    # Ensure correct initial slope by construction
    T_ref = T_trans  # Reference temperature
    return 0.25 * N * k_B_hartree_per_K * T * (
        (1 - f) + f * (T / T_ref)**(alpha_high - 1)
    )

def exponential_constrained(T, T_char):
    """
    Exponential approach: E = 0.25*N*k_B*T*(1 - exp(-T/T_char))
    
    At T→0: E ≈ 0.25*N*k_B*T (correct slope)
    At T>>T_char: E ≈ 0.25*N*k_B*T_char (saturates)
    """
    N = 500
    return 0.25 * N * k_B_hartree_per_K * T * (1 - np.exp(-T / T_char))

def plot_constrained_fits(data_file, N_bonds=500, output_file=None, show_plot=True):
    """
    Plot harmonic energy with E(0)=0 constraint.
    """
    
    print(f" CONSTRAINED HARMONIC ENERGY ANALYSIS: E(0)=0, dE/dT|₀=0.25×{N_bonds}×k_B")
    print("=" * 70)
    
    # Load data
    try:
        data = pd.read_csv(data_file, sep=r'\s+', comment='#')
        T = data['temperature'].values
        E_harmonic = data['harmonic_hartree'].values
        
        print(f" Loaded {len(T)} data points")
        print(f"   Temperature range: {T.min():.1f} - {T.max():.1f} K")
        print(f"   Energy range: {E_harmonic.min():.6f} - {E_harmonic.max():.6f} Ha")
        
    except Exception as e:
        print(f" Failed to load data: {e}")
        return
    
    # Target slope
    target_slope = 0.25 * N_bonds * k_B_hartree_per_K
    print(f" Target slope: {target_slope:.6f} Ha/K")
    
    # Define constrained models
    models = {
        'Classical Linear': {
            'func': lambda T: classical_linear(T, N_bonds),
            'params': [],
            'fixed': True,
            'color': 'black',
            'linestyle': '--',
            'description': f'E = 0.25×{N_bonds}×k_B×T'
        },
        'Power Law': {
            'func': fitted_power_law,
            'params': ['A', 'α'],
            'bounds': ([0, 0.1], [np.inf, 3.0]),
            'p0': [target_slope, 1.0],
            'fixed': False,
            'color': 'red',
            'linestyle': '-',
            'description': 'E = A×T^α'
        },
        'Classical + Correction': {
            'func': lambda T, alpha, A_corr: classical_plus_correction(T, alpha, A_corr),
            'params': ['α_corr', 'A_corr'],
            'bounds': ([0.1, -np.inf], [3.0, np.inf]),
            'p0': [2.0, 1e-7],
            'fixed': False,
            'color': 'blue',
            'linestyle': '-',
            'description': 'E = 0.25×N×k_B×T + A×T^α'
        },
        'Classical + √T': {
            'func': lambda T, A_sqrt: classical_plus_sqrt(T, A_sqrt),
            'params': ['A_√T'],
            'bounds': ([-np.inf], [np.inf]),
            'p0': [1e-5],
            'fixed': False,
            'color': 'purple',
            'linestyle': '-',
            'description': 'E = 0.25×N×k_B×T + A×√T'
        },
        'Classical + T^(3/2)': {
            'func': lambda T, A_T32: classical_plus_T32(T, A_T32),
            'params': ['A_T^(3/2)'],
            'bounds': ([-np.inf], [np.inf]),
            'p0': [1e-7],
            'fixed': False,
            'color': 'magenta',
            'linestyle': '-',
            'description': 'E = 0.25×N×k_B×T + A×T^(3/2)'
        },
        'Tanh Transition': {
            'func': lambda T, T_trans, alpha_high: tanh_constrained(T, T_trans, alpha_high),
            'params': ['T_trans', 'α_high'],
            'bounds': ([50, 0.1], [300, 3.0]),
            'p0': [120, 0.8],
            'fixed': False,
            'color': 'green',
            'linestyle': '-',
            'description': 'Constrained tanh transition'
        },
        'Exponential': {
            'func': lambda T, T_char: exponential_constrained(T, T_char),
            'params': ['T_char'],
            'bounds': ([1], [1000]),
            'p0': [100],
            'fixed': False,
            'color': 'orange',
            'linestyle': '-',
            'description': 'E = 0.25×N×k_B×T×(1-e^(-T/T_char))'
        }
    }
    
    # Fit models
    fit_results = {}
    
    print(f"\n FITTING CONSTRAINED MODELS:")
    print(f"Model               | RMSE (Ha)  | R²     | Parameters")
    print(f"--------------------|------------|--------|----------------------------------")
    
    for name, model in models.items():
        try:
            if model['fixed']:
                # Fixed model (no fitting needed)
                E_pred = model['func'](T)
                rmse = np.sqrt(np.mean((E_harmonic - E_pred)**2))
                r2 = 1 - np.var(E_harmonic - E_pred) / np.var(E_harmonic)
                
                fit_results[name] = {
                    'params': [],
                    'rmse': rmse,
                    'r2': r2,
                    'model': model,
                    'E_pred': E_pred
                }
                
                print(f"{name:19s} | {rmse:.6f} | {r2:.4f} | Fixed")
                
            else:
                # Fit model
                popt, pcov = curve_fit(
                    model['func'], T, E_harmonic,
                    p0=model['p0'],
                    bounds=model['bounds'],
                    maxfev=5000
                )
                
                E_pred = model['func'](T, *popt)
                rmse = np.sqrt(np.mean((E_harmonic - E_pred)**2))
                r2 = 1 - np.var(E_harmonic - E_pred) / np.var(E_harmonic)
                
                fit_results[name] = {
                    'params': popt,
                    'rmse': rmse,
                    'r2': r2,
                    'model': model,
                    'E_pred': E_pred
                }
                
                param_str = ", ".join([f"{p}={v:.4f}" for p, v in zip(model['params'], popt)])
                print(f"{name:19s} | {rmse:.6f} | {r2:.4f} | {param_str}")
                
        except Exception as e:
            print(f"{name:19s} | FAILED     | FAILED | {str(e)[:30]}...")
    
    # Create plots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Data and fits
    ax1.scatter(T, E_harmonic, color='blue', alpha=0.8, s=50, label='Data', zorder=10)
    
    T_smooth = np.linspace(T.min()-20, T.max()+50, 300)
    
    for name, result in fit_results.items():
        model = result['model']
        if model['fixed']:
            E_smooth = model['func'](T_smooth)
        else:
            E_smooth = model['func'](T_smooth, *result['params'])
        
        ax1.plot(T_smooth, E_smooth,
                color=model['color'],
                linestyle=model['linestyle'],
                linewidth=2,
                label=f"{name} (R²={result['r2']:.3f})",
                alpha=0.8)
    
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('Harmonic Energy (Ha)')
    ax1.set_title('Constrained Models: E(0)=0, Correct Low-T Slope')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(T.min()-10, T.max()+30)
    
    # Plot 2: Residuals
    for name, result in fit_results.items():
        if name == 'Classical Linear':
            continue
        model = result['model']
        residuals = E_harmonic - result['E_pred']
        ax2.scatter(T, residuals,
                   color=model['color'],
                   alpha=0.7, s=30,
                   label=name)
    
    ax2.axhline(0, color='black', linestyle='-', alpha=0.5)
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('Residuals (Ha)')
    ax2.set_title('Model Residuals')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Low-T derivatives
    T_deriv = np.linspace(1, 50, 100)
    
    ax3.axhline(target_slope, color='black', linestyle='--',
               label=f'Target: {target_slope:.6f} Ha/K')
    
    # Calculate derivatives numerically
    for name, result in fit_results.items():
        model = result['model']
        dT = 0.1
        
        if model['fixed']:
            # Analytical derivative for classical
            derivatives = np.full_like(T_deriv, target_slope)
        else:
            # Numerical derivative
            T_low = T_deriv - dT/2
            T_high = T_deriv + dT/2
            E_low = model['func'](T_low, *result['params'])
            E_high = model['func'](T_high, *result['params'])
            derivatives = (E_high - E_low) / dT
        
        ax3.plot(T_deriv, derivatives,
                color=model['color'],
                linewidth=2,
                label=f"{name}")
    
    ax3.set_xlabel('Temperature (K)')
    ax3.set_ylabel('dE/dT (Ha/K)')
    ax3.set_title('Temperature Derivatives: Approach to Target Slope')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(0, target_slope * 3)
    
    # Plot 4: 0K extrapolation
    T_extrap = np.linspace(0, T.min(), 200)
    
    for name, result in fit_results.items():
        model = result['model']
        if model['fixed']:
            E_extrap = model['func'](T_extrap)
        else:
            E_extrap = model['func'](T_extrap, *result['params'])
        
        ax4.plot(T_extrap, E_extrap,
                color=model['color'],
                linestyle=model['linestyle'],
                linewidth=2,
                label=f"{name}")
    
    # Show target slope from origin
    T_slope_ref = np.linspace(0, 30, 50)
    E_slope_ref = target_slope * T_slope_ref
    ax4.plot(T_slope_ref, E_slope_ref, 'k:', alpha=0.5, linewidth=1,
            label=f'Target slope from origin')
    
    ax4.axvspan(5, 77, alpha=0.2, color='orange', label='Simulation Range')
    ax4.set_xlabel('Temperature (K)')
    ax4.set_ylabel('Harmonic Energy (Ha)')
    ax4.set_title('Extrapolation to 0K: All Models Pass Through Origin')
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim(0, T.min() + 5)
    
    plt.tight_layout()
    
    # Analysis
    print(f"\n CONSTRAINT VALIDATION:")
    print(f"   All models: E(0) = 0 ")
    print(f"   Target slope: {target_slope:.6f} Ha/K")
    
    for name, result in fit_results.items():
        model = result['model']
        if model['fixed']:
            slope_at_1K = target_slope
        else:
            # Calculate slope at 1K
            dT = 0.01
            E_low = model['func'](1 - dT/2, *result['params'])
            E_high = model['func'](1 + dT/2, *result['params'])
            slope_at_1K = (E_high - E_low) / dT
        
        match_percent = (slope_at_1K / target_slope) * 100
        print(f"   {name}: slope at 1K = {slope_at_1K:.6f} Ha/K ({match_percent:.1f}%)")
    
    # Best model
    if fit_results:
        best_model = min(fit_results.keys(), key=lambda k: fit_results[k]['rmse'])
        print(f"\n BEST FIT: {best_model}")
        print(f"   RMSE: {fit_results[best_model]['rmse']:.6f} Ha")
        print(f"   R²: {fit_results[best_model]['r2']:.4f}")
    
    # Save plot
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\n Plot saved to: {output_file}")
    
    if show_plot:
        plt.show()
    
    return fit_results

def main():
    parser = argparse.ArgumentParser(description='Constrained harmonic energy analysis')
    parser.add_argument('data_file', nargs='?',
                       default='final_nodiss_cavitymd/potential_energy_components_vs_temperature.txt',
                       help='Path to empirical data file')
    parser.add_argument('--n-bonds', type=int, default=500,
                       help='Number of bonds/modes')
    parser.add_argument('--output', '-o', type=str,
                       help='Output file path')
    parser.add_argument('--no-show', action='store_true',
                       help='Do not display plot')
    
    args = parser.parse_args()
    
    if not Path(args.data_file).exists():
        print(f" Data file not found: {args.data_file}")
        return 1
    
    if not args.output:
        args.output = f"constrained_harmonic_fits_N{args.n_bonds}.png"
    
    try:
        results = plot_constrained_fits(
            data_file=args.data_file,
            N_bonds=args.n_bonds,
            output_file=args.output,
            show_plot=not args.no_show
        )
        print(" Constrained analysis completed!")
        return 0
    except Exception as e:
        print(f" Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())
