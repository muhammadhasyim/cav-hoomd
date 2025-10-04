#!/usr/bin/env python3
"""
Plot Harmonic Energy vs Temperature with Different Correction Models

This script plots harmonic energy data and compares various models for 
deviations from the classical linear temperature dependence.
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
    Exact classical linear: E = 0.25 * N * k_B * T
    
    Constraints: E(0) = 0 and dE/dT = 0.25 * N * k_B
    """
    return 0.25 * N * k_B_hartree_per_K * T

def power_law_correction(T, alpha, N=500):
    """
    Power law: E = A * T^α where A chosen to give correct slope at T→0
    
    For dE/dT|_{T→0} = 0.25 * N * k_B, we need:
    A * α = 0.25 * N * k_B, so A = 0.25 * N * k_B / α
    
    But this only works for α = 1. For α ≠ 1, we fit A directly.
    """
    if alpha == 1.0:
        # Exact classical case
        return 0.25 * N * k_B_hartree_per_K * T
    else:
        # For α ≠ 1, we need to fit A to match the data
        # The constraint dE/dT|_{T→0} = 0.25*N*k_B only works for α=1
        # For α < 1, dE/dT|_{T→0} = 0, for α > 1, dE/dT|_{T→0} = ∞
        # So we'll fit A without the constraint for α ≠ 1
        return (0.25 * N * k_B_hartree_per_K) * (T**alpha)

def exponential_correction(T, E0, N, T_char):
    """
    Exponential correction for temperature-dependent effective modes:
    E = E₀ + 0.25 * N * k_B * T * (1 - exp(-T/T_char))
    
    At low T: E ≈ E₀ + 0.25 * N * k_B * T²/T_char (quadratic)
    At high T: E ≈ E₀ + 0.25 * N * k_B * T (linear, correct slope)
    """
    return E0 + 0.25 * N * k_B_hartree_per_K * T * (1 - np.exp(-T/T_char))

def coupling_correction(T, E0, N, T_coupling):
    """
    Coupling to environment correction:
    E = E₀ + 0.25 * N * k_B * T / (1 + T/T_coupling)
    
    At low T: E ≈ E₀ + 0.25 * N * k_B * T (linear, correct slope)
    At high T: E ≈ E₀ + 0.25 * N * k_B * T_coupling (saturates)
    """
    return E0 + 0.25 * N * k_B_hartree_per_K * T / (1 + T/T_coupling)

def logarithmic_correction(T, E0, N, T_scale):
    """
    Logarithmic correction for frustrated modes:
    E = E₀ + 0.25 * N * k_B * T_scale * ln(1 + T/T_scale)
    
    At low T: E ≈ E₀ + 0.25 * N * k_B * T (linear, correct slope)
    At high T: E ≈ E₀ + 0.25 * N * k_B * T_scale * ln(T) (logarithmic growth)
    """
    return E0 + 0.25 * N * k_B_hartree_per_K * T_scale * np.log(1 + T/T_scale)

def combined_correction(T, E0, N, alpha, A_extra):
    """
    Combined model: Linear + power law correction
    E = E₀ + 0.25 * N * k_B * T + A_extra * T^α
    """
    return E0 + 0.25 * N * k_B_hartree_per_K * T + A_extra * (T**alpha)

def low_T_classical(T, E0, N, T_cross):
    """
    Low-T classical model: Approaches linear at low T, power law at high T
    E = E₀ + N * k_B * T * (T/(T + T_cross))
    
    At low T (T << T_cross): E ≈ E₀ + N * k_B * T²/T_cross (quadratic → can be made linear)
    At high T (T >> T_cross): E ≈ E₀ + N * k_B * T (linear, classical)
    """
    return E0 + N * k_B_hartree_per_K * T * (T / (T + T_cross))

def piecewise_linear(T, E0, N, T_break, slope_high):
    """
    Piecewise model: Linear at low T, different slope at high T
    Low T:  E = E₀ + 0.25 * N * k_B * T
    High T: E = E₀ + 0.25 * N * k_B * T_break + slope_high * (T - T_break)
    """
    result = np.zeros_like(T)
    low_mask = T <= T_break
    high_mask = T > T_break
    
    # Low temperature: classical linear with correct slope
    result[low_mask] = E0 + 0.25 * N * k_B_hartree_per_K * T[low_mask]
    
    # High temperature: different slope
    E_at_break = E0 + 0.25 * N * k_B_hartree_per_K * T_break
    result[high_mask] = E_at_break + slope_high * (T[high_mask] - T_break)
    
    return result

def tanh_transition(T, E0, N, T_trans, alpha_high):
    """
    Smooth transition from classical linear to power law behavior
    Uses tanh function to smoothly interpolate between regimes
    
    Low T: E ≈ E₀ + 0.25 * N * k_B * T (classical linear)  
    High T: E ≈ E₀ + 0.25 * N * k_B * T^α_high (power law)
    """
    # Smooth interpolation factor
    f = 0.5 * (1 + np.tanh((T - T_trans) / (0.2 * T_trans)))
    
    # Linear component (dominant at low T)
    E_linear = E0 + 0.25 * N * k_B_hartree_per_K * T
    
    # Power law component (dominant at high T)  
    E_power = E0 + 0.25 * N * k_B_hartree_per_K * (T**alpha_high)
    
    # Smooth interpolation
    return (1 - f) * E_linear + f * E_power

def crossover_model(T, E0, N_eff, T_char, alpha):
    """
    Crossover model: Classical at low T, sublinear at high T
    E = E₀ + 0.25 * N_eff * k_B * T_char * (T/T_char)^α * tanh(T/T_char)
    
    At low T: E ≈ E₀ + 0.25 * N_eff * k_B * T (classical linear)
    At high T: E ≈ E₀ + 0.25 * N_eff * k_B * T_char * (T/T_char)^α (power law)
    """
    x = T / T_char
    return E0 + 0.25 * N_eff * k_B_hartree_per_K * T_char * (x**alpha) * np.tanh(x)

def plot_harmonic_energy_fits(data_file, N_bonds=500, output_file=None, show_plot=True):
    """
    Plot harmonic energy vs temperature with different correction models.
    
    Parameters
    ----------
    data_file : str
        Path to empirical energy-temperature data file
    N_bonds : int
        Number of harmonic bonds/modes (default: 500)
    output_file : str, optional
        Output file for saving the plot
    show_plot : bool, optional
        Whether to display the plot (default: True)
    """
    
    print(f" HARMONIC ENERGY ANALYSIS: {N_bonds} bonds/modes")
    print("=" * 60)
    
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
    
    # Define models to fit
    models = {
        'Classical Linear': {
            'func': lambda T, E0: classical_linear(T, E0, N_bonds),
            'params': ['E0'],
            'bounds': ([-np.inf], [np.inf]),
            'p0': [E_harmonic.min()],
            'color': 'black',
            'linestyle': '--',
            'description': f'E = E₀ + {N_bonds}×k_B×T'
        },
        'Power Law': {
            'func': lambda T, E0, alpha: power_law_correction(T, E0, N_bonds, alpha),
            'params': ['E0', 'α'],
            'bounds': ([-np.inf, 0.1], [np.inf, 3.0]),
            'p0': [E_harmonic.min(), 1.0],
            'color': 'red',
            'linestyle': '-',
            'description': f'E = E₀ + {N_bonds}×k_B×T^α'
        },
        'Exponential': {
            'func': lambda T, E0, T_char: exponential_correction(T, E0, N_bonds, T_char),
            'params': ['E0', 'T_char'],
            'bounds': ([-np.inf, 1], [np.inf, 1000]),
            'p0': [E_harmonic.min(), 100],
            'color': 'blue',
            'linestyle': '-',
            'description': f'E = E₀ + {N_bonds}×k_B×T×(1-e^(-T/T_char))'
        },
        'Coupling': {
            'func': lambda T, E0, T_coupling: coupling_correction(T, E0, N_bonds, T_coupling),
            'params': ['E0', 'T_coupling'],
            'bounds': ([-np.inf, 1], [np.inf, 10000]),
            'p0': [E_harmonic.min(), 300],
            'color': 'green',
            'linestyle': '-',
            'description': f'E = E₀ + {N_bonds}×k_B×T/(1+T/T_coupling)'
        },
        'Logarithmic': {
            'func': lambda T, E0, T_scale: logarithmic_correction(T, E0, N_bonds, T_scale),
            'params': ['E0', 'T_scale'],
            'bounds': ([-np.inf, 1], [np.inf, 1000]),
            'p0': [E_harmonic.min(), 100],
            'color': 'orange',
            'linestyle': '-',
            'description': f'E = E₀ + {N_bonds}×k_B×T_scale×ln(1+T/T_scale)'
        },
        'Combined': {
            'func': lambda T, E0, alpha, A_extra: combined_correction(T, E0, N_bonds, alpha, A_extra),
            'params': ['E0', 'α', 'A_extra'],
            'bounds': ([-np.inf, 0.1, -np.inf], [np.inf, 3.0, np.inf]),
            'p0': [E_harmonic.min(), 0.8, 1e-6],
            'color': 'purple',
            'linestyle': '-',
            'description': f'E = E₀ + {N_bonds}×k_B×T + A×T^α'
        },
        'Tanh Transition': {
            'func': lambda T, E0, T_trans, alpha_high: tanh_transition(T, E0, N_bonds, T_trans, alpha_high),
            'params': ['E0', 'T_trans', 'α_high'],
            'bounds': ([-np.inf, 10, 0.1], [np.inf, 500, 2.0]),
            'p0': [E_harmonic.min(), 100, 0.8],
            'color': 'cyan',
            'linestyle': '-',
            'description': f'Linear→Power transition at T_trans'
        },
        'Crossover': {
            'func': lambda T, E0, N_eff, T_char, alpha: crossover_model(T, E0, N_eff, T_char, alpha),
            'params': ['E0', 'N_eff', 'T_char', 'α'],
            'bounds': ([-np.inf, 10, 1, 0.1], [np.inf, N_bonds*2, 1000, 2.0]),
            'p0': [E_harmonic.min(), N_bonds//2, 100, 0.8],
            'color': 'magenta',
            'linestyle': '-',
            'description': f'Classical→sublinear crossover'
        },
        'Piecewise': {
            'func': lambda T, E0, T_break, slope_high: piecewise_linear(T, E0, N_bonds, T_break, slope_high),
            'params': ['E0', 'T_break', 'slope_high'],
            'bounds': ([-np.inf, 50, 0], [np.inf, 300, 0.01]),
            'p0': [E_harmonic.min(), 120, 0.0001],
            'color': 'brown',
            'linestyle': '-',
            'description': f'Piecewise: classical→different slope'
        }
    }
    
    # Fit all models
    fit_results = {}
    
    print(f"\n FITTING MODELS:")
    print(f"Model           | RMSE (Ha)  | R²     | Parameters")
    print(f"----------------|------------|--------|----------------------------------")
    
    for name, model in models.items():
        try:
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
                'model': model
            }
            
            # Format parameters for display
            param_str = ", ".join([f"{p}={v:.4f}" for p, v in zip(model['params'], popt)])
            print(f"{name:15s} | {rmse:.6f} | {r2:.4f} | {param_str}")
            
        except Exception as e:
            print(f"{name:15s} | FAILED     | FAILED | {str(e)[:20]}...")
            continue
    
    # Create the plot
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Data and all model fits
    ax1.scatter(T, E_harmonic, color='blue', alpha=0.8, s=50, label='Data', zorder=10)
    
    T_smooth = np.linspace(T.min()-20, T.max()+50, 300)
    
    for name, result in fit_results.items():
        model = result['model']
        E_smooth = model['func'](T_smooth, *result['params'])
        ax1.plot(T_smooth, E_smooth, 
                color=model['color'], 
                linestyle=model['linestyle'],
                linewidth=2,
                label=f"{name} (R²={result['r2']:.3f})",
                alpha=0.8)
    
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('Harmonic Energy (Ha)')
    ax1.set_title('Harmonic Energy vs Temperature: Model Comparison')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(T.min()-10, T.max()+30)
    
    # Plot 2: Residuals
    for name, result in fit_results.items():
        if name == 'Classical Linear':
            continue  # Skip classical for residuals plot
        model = result['model']
        E_pred = model['func'](T, *result['params'])
        residuals = E_harmonic - E_pred
        ax2.scatter(T, residuals, 
                   color=model['color'], 
                   alpha=0.7, s=30,
                   label=f"{name}")
    
    ax2.axhline(0, color='black', linestyle='-', alpha=0.5)
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('Residuals (Ha)')
    ax2.set_title('Model Residuals')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Effective temperature dependence (derivatives)
    ax3.set_xlabel('Temperature (K)')
    ax3.set_ylabel('dE/dT (Ha/K)')
    ax3.set_title('Temperature Dependence: dE/dT')
    
    # Classical expectation
    classical_slope = N_bonds * k_B_hartree_per_K
    ax3.axhline(classical_slope, color='black', linestyle='--', 
               label=f'Classical: {classical_slope:.6f} Ha/K')
    
    # Plot derivatives for fitted models
    T_deriv = np.linspace(T.min(), T.max()+100, 200)
    
    # Power law derivative
    if 'Power Law' in fit_results:
        popt = fit_results['Power Law']['params']
        E0, alpha = popt
        derivative = N_bonds * k_B_hartree_per_K * alpha * (T_deriv**(alpha-1))
        ax3.plot(T_deriv, derivative, color='red', linewidth=2, 
                label=f'Power Law (α={alpha:.3f})')
    
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(T.min(), T.max()+50)
    
    # Plot 4: Low-temperature extrapolation (0K to data minimum)
    T_low = np.linspace(0, T.min(), 200)
    
    # Plot models that should show good low-T behavior
    priority_models = ['Classical Linear', 'Power Law', 'Tanh Transition', 'Crossover', 'Combined']
    
    for name, result in fit_results.items():
        if name in priority_models:
            model = result['model']
            E_low = model['func'](T_low, *result['params'])
            ax4.plot(T_low, E_low, 
                    color=model['color'],
                    linestyle=model['linestyle'],
                    linewidth=2,
                    label=f"{name} (R²={result['r2']:.3f})")
    
    # Add reference lines
    ax4.axvspan(5, 77, alpha=0.2, color='orange', label='Simulation Range')
    ax4.axvspan(0, 20, alpha=0.1, color='lightblue', label='Deep Extrapolation')
    
    # Show classical slope at very low T
    classical_slope = 0.25 * N_bonds * k_B_hartree_per_K
    if 'Classical Linear' in fit_results:
        E0_classical = fit_results['Classical Linear']['params'][0]
        T_slope_demo = np.linspace(0, 30, 50)
        E_slope_demo = E0_classical + classical_slope * T_slope_demo
        ax4.plot(T_slope_demo, E_slope_demo, 'k:', alpha=0.5, linewidth=1,
                label=f'Target slope: {classical_slope:.6f} Ha/K')
    
    ax4.set_xlabel('Temperature (K)')
    ax4.set_ylabel('Harmonic Energy (Ha)')
    ax4.set_title('Low-Temperature Extrapolation: 0K → Data Range')
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim(0, T.min() + 5)
    
    plt.tight_layout()
    
    # Physical analysis
    print(f"\n PHYSICAL ANALYSIS:")
    
    if 'Classical Linear' in fit_results and 'Power Law' in fit_results:
        classical_params = fit_results['Classical Linear']['params']
        power_params = fit_results['Power Law']['params']
        
        E0_classical = classical_params[0]
        E0_power, alpha = power_params
        
        print(f"   Classical E₀: {E0_classical:.6f} Ha")
        print(f"   Power law E₀: {E0_power:.6f} Ha, α = {alpha:.3f}")
        
        # Check what this implies about the effective number of modes
        fitted_slope_at_200K = 0.25 * N_bonds * k_B_hartree_per_K * alpha * (200**(alpha-1))
        effective_N_at_200K = fitted_slope_at_200K / (0.25 * k_B_hartree_per_K)
        
        print(f"   Effective modes at 200K: {effective_N_at_200K:.0f} (vs {N_bonds} assumed)")
        print(f"   Mode utilization: {effective_N_at_200K/N_bonds:.1%}")
        
        if alpha < 0.9:
            print(f"     Strong sublinear behavior suggests significant coupling/constraints")
        elif alpha > 1.1:
            print(f"     Superlinear behavior suggests anharmonic effects")
        else:
            print(f"    Nearly classical behavior")
    
    # Analyze low-temperature behavior for transition models
    print(f"\n  LOW-TEMPERATURE BEHAVIOR ANALYSIS:")
    
    T_test_low = np.array([1, 5, 10, 20])
    classical_slope = 0.25 * N_bonds * k_B_hartree_per_K
    
    print(f"   Target low-T slope: {classical_slope:.6f} Ha/K (0.25 × {N_bonds} × k_B)")
    print(f"   Model slopes at low T:")
    print(f"   T(K) | Classical | Power | Tanh | Crossover")
    print(f"   -----|-----------|-------|------|----------")
    
    for T_test in T_test_low:
        slopes = {}
        
        # Classical (constant)
        slopes['Classical'] = classical_slope
        
        # Power law derivative
        if 'Power Law' in fit_results:
            popt = fit_results['Power Law']['params']
            alpha = popt[1] if len(popt) > 1 else 1.0
            slopes['Power'] = 0.25 * N_bonds * k_B_hartree_per_K * alpha * (T_test**(alpha-1))
        
        # Tanh transition - calculate numerical derivative
        if 'Tanh Transition' in fit_results:
            popt = fit_results['Tanh Transition']['params']
            dT = 0.1
            E1 = tanh_transition(T_test - dT/2, popt[0], N_bonds, popt[1], popt[2])
            E2 = tanh_transition(T_test + dT/2, popt[0], N_bonds, popt[1], popt[2])
            slopes['Tanh'] = (E2 - E1) / dT
            
        # Crossover model - calculate numerical derivative  
        if 'Crossover' in fit_results:
            popt = fit_results['Crossover']['params']
            dT = 0.1
            E1 = crossover_model(T_test - dT/2, popt[0], popt[1], popt[2], popt[3])
            E2 = crossover_model(T_test + dT/2, popt[0], popt[1], popt[2], popt[3])
            slopes['Crossover'] = (E2 - E1) / dT
        
        # Print slopes
        slope_strs = []
        for model in ['Classical', 'Power', 'Tanh', 'Crossover']:
            if model in slopes:
                ratio = slopes[model] / classical_slope
                slope_strs.append(f"{slopes[model]:.6f}")
            else:
                slope_strs.append("N/A")
        
        print(f"   {T_test:4.0f} | {slope_strs[0]} | {slope_strs[1]} | {slope_strs[2]} | {slope_strs[3]}")
    
    # Check which models approach classical behavior at low T
    print(f"\n    Models approaching classical slope at low T:")
    for model_name in ['Tanh Transition', 'Crossover']:
        if model_name in fit_results:
            # Calculate slope at T=1K
            if model_name == 'Tanh Transition':
                popt = fit_results[model_name]['params']
                dT = 0.01
                E1 = tanh_transition(1 - dT/2, popt[0], N_bonds, popt[1], popt[2])
                E2 = tanh_transition(1 + dT/2, popt[0], N_bonds, popt[1], popt[2])
                slope_1K = (E2 - E1) / dT
            elif model_name == 'Crossover':
                popt = fit_results[model_name]['params']
                dT = 0.01
                E1 = crossover_model(1 - dT/2, popt[0], popt[1], popt[2], popt[3])
                E2 = crossover_model(1 + dT/2, popt[0], popt[1], popt[2], popt[3])
                slope_1K = (E2 - E1) / dT
            
            ratio_1K = slope_1K / classical_slope
            print(f"     {model_name}: {ratio_1K:.1%} of classical at 1K")
            
            if ratio_1K > 0.8:
                print(f"        Good low-T classical behavior")
            elif ratio_1K > 0.5:
                print(f"         Moderate approach to classical")
            else:
                print(f"        Poor low-T classical behavior")
    
    # Best model recommendation
    if fit_results:
        best_model = min(fit_results.keys(), key=lambda k: fit_results[k]['rmse'])
        print(f"\n BEST FIT: {best_model}")
        print(f"   RMSE: {fit_results[best_model]['rmse']:.6f} Ha")
        print(f"   R²: {fit_results[best_model]['r2']:.4f}")
        print(f"   Description: {fit_results[best_model]['model']['description']}")
    
    # Save plot
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\n Plot saved to: {output_file}")
    
    if show_plot:
        plt.show()
    
    return fit_results

def main():
    parser = argparse.ArgumentParser(
        description='Plot harmonic energy vs temperature with correction models',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('data_file', nargs='?', 
                       default='final_nodiss_cavitymd/potential_energy_components_vs_temperature.txt',
                       help='Path to empirical energy-temperature data file')
    parser.add_argument('--n-bonds', type=int, default=500,
                       help='Number of harmonic bonds/modes (default: 500)')
    parser.add_argument('--output', '-o', type=str,
                       help='Output file path for saving the plot')
    parser.add_argument('--no-show', action='store_true',
                       help='Do not display the plot')
    
    args = parser.parse_args()
    
    # Check if file exists
    if not Path(args.data_file).exists():
        print(f" Data file not found: {args.data_file}")
        return 1
    
    # Generate output filename if not provided
    if not args.output:
        args.output = f"harmonic_energy_fits_N{args.n_bonds}.png"
    
    # Run the analysis
    try:
        results = plot_harmonic_energy_fits(
            data_file=args.data_file,
            N_bonds=args.n_bonds,
            output_file=args.output,
            show_plot=not args.no_show
        )
        print(" Analysis completed successfully!")
        return 0
        
    except Exception as e:
        print(f" Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())
