#!/usr/bin/env python3
"""
Run Material Time Analysis with CORRECTED Reconstruction Method

This script uses the corrected material time reconstruction from material_time_correct.py
and generates the stretched exponential fits and unified collapse plots.

Author: Material Time Analysis (Corrected)
Date: November 2, 2025
"""

import os
import sys
import numpy as np
from pathlib import Path

REPRO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPRO_ROOT / "shared"))

# Import the corrected reconstructor
from material_time_correct import MaterialTimeReconstructor

# Import the old analyzer for its plotting methods
from material_time_analysis import MaterialTimeAnalyzer


class CorrectedMaterialTimeAnalyzer(MaterialTimeAnalyzer):
    """
    Material Time Analyzer using the CORRECTED reconstruction method.
    
    Inherits all plotting and analysis methods from MaterialTimeAnalyzer,
    but uses the corrected least squares reconstruction for material time.
    """
    
    def __init__(self, criterion_value=0.1, alpha_smooth=1e-10, n_grid_points=500, 
                 sg_window_ps=0.5, verbose=True):
        """
        Initialize with corrected reconstructor parameters.
        
        Parameters:
        -----------
        criterion_value : float
            Criterion for material time unit (default 0.1)
        alpha_smooth : float
            Regularization parameter for smoothness
        n_grid_points : int
            Number of grid points for reconstruction
        sg_window_ps : float
            Savitzky-Golay filter window (for compatibility)
        verbose : bool
            Verbose output
        """
        super().__init__(criterion_value=criterion_value, sg_window_ps=sg_window_ps, verbose=verbose)
        
        # Initialize the corrected reconstructor
        self.reconstructor = MaterialTimeReconstructor(
            criterion_value=criterion_value,
            alpha_smooth=alpha_smooth,
            n_grid_points=n_grid_points,
            verbose=False  # Reduce verbosity for batch processing
        )
    
    def calculate_material_time(self, normalized_fkt):
        """
        Calculate material time using the CORRECTED least squares method.
        
        This overrides the old (incorrect) sequential interpolation method
        with the new simultaneous constraint solving approach.
        
        Parameters:
        -----------
        normalized_fkt : dict
            Normalized F(k,t) data for all refs
            
        Returns:
        --------
        time_grid : array
            Laboratory time grid
        xi_grid : array
            Material time ξ(t) on the grid
        """
        if self.verbose:
            print("\n" + "="*70)
            print("Using CORRECTED material time reconstruction")
            print("Method: Simultaneous least squares with regularization")
            print(f"  α = {self.reconstructor.alpha_smooth}")
            print(f"  Grid points = {self.reconstructor.n_grid_points}")
            print("="*70)
        
        # Convert normalized_fkt to the format expected by the reconstructor
        # normalized_fkt format: ref_num -> (lag_time, fkt)
        
        # Calculate waiting times (same as before)
        waiting_times = {}
        for ref_num in normalized_fkt.keys():
            if ref_num == 0 or ref_num == 14:
                continue
            waiting_times[ref_num] = (ref_num - 1) * 200
        
        # Find crossings
        crossings = self.reconstructor.find_crossings(normalized_fkt, waiting_times)
        
        if len(crossings) == 0:
            if self.verbose:
                print("WARNING: No crossings found! Using linear fallback")
            time_grid = np.linspace(0, 2800, 1000)
            xi_grid = time_grid / np.max(time_grid) * 2.0
            return time_grid, xi_grid
        
        # Build and solve constraint system
        time_points, A, b = self.reconstructor.build_constraint_system(crossings)
        xi_values = self.reconstructor.solve_least_squares(time_points, A, b)
        
        # Interpolate to fine grid
        time_grid, xi_grid = self.reconstructor.interpolate_to_fine_grid(time_points, xi_values)
        
        if self.verbose:
            print(f"\nCorrected reconstruction complete:")
            print(f"  ξ range: [{xi_grid.min():.3f}, {xi_grid.max():.3f}]")
            print(f"  Time range: [0, {time_grid[-1]:.1f}] ps")
        
        return time_grid, xi_grid


def run_corrected_analysis(base_dir, profile=None):
    """
    Run the complete material time analysis with corrected reconstruction.

    Generates:
    - stretched_exponential_fits.png
    - unified_material_time_collapse.png
    - overlapped_stretched_exponential_fits.png

    Parameters
    ----------
    base_dir : str or Path
        Output directory for plots and derived tables.
    profile : DatasetProfile, optional
        When set, coupling directories are taken from the staged layout.
    """
    print("\n" + "="*80)
    print("MATERIAL TIME ANALYSIS WITH CORRECTED RECONSTRUCTION")
    print("="*80)
    print("\nKey improvements:")
    print("  • Simultaneous constraint solving (not sequential interpolation)")
    print("  • Regularized least squares for smoothness")
    print("  • Dense regular grid (500 points) for smooth curves")
    print("="*80)
    
    base_dir = Path(base_dir)

    if profile is not None:
        coupling_dirs = [
            str(profile.staged_root / entry.staged_dir_name)
            for entry in profile.couplings
        ]
    else:
        coupling_dirs = [
            os.path.join(base_dir, "cavity_coupling_0epos00_switch_200.0ps"),
            os.path.join(base_dir, "cavity_coupling_3eneg04_switch_200.0ps"),
            os.path.join(base_dir, "cavity_coupling_5eneg04_switch_200.0ps"),
            os.path.join(base_dir, "cavity_coupling_7eneg04_switch_200.0ps"),
            os.path.join(base_dir, "cavity_coupling_1eneg03_switch_200.0ps"),
        ]
    
    # Initialize the corrected analyzer
    analyzer = CorrectedMaterialTimeAnalyzer(
        criterion_value=0.1,
        alpha_smooth=0.01,  # Very small regularization
        n_grid_points=50,
        verbose=True
    )
    
    # Process each coupling strength
    successful_analyses = []
    
    for coupling_dir in coupling_dirs:
        coupling_name = os.path.basename(coupling_dir)
        
        if not os.path.exists(coupling_dir):
            print(f"\nWarning: {coupling_dir} does not exist, skipping")
            continue
        
        print(f"\n{'='*70}")
        print(f"Processing: {coupling_name}")
        print(f"{'='*70}")
        
        try:
            # Load and normalize F(k,t) data
            fkt_files = analyzer.load_fkt_data(coupling_dir)
            waiting_times = analyzer.calculate_waiting_times(fkt_files)
            normalized_fkt = analyzer.normalize_fkt_data(fkt_files)
            
            # Calculate material time using CORRECTED method
            time_grid, xi_grid = analyzer.calculate_material_time(normalized_fkt)
            
            # Collapse data
            collapsed_data = analyzer.collapse_data(
                (time_grid, xi_grid), 
                normalized_fkt, 
                waiting_times
            )
            
            # Store results
            analyzer.collapsed_data[coupling_dir] = collapsed_data
            successful_analyses.append(coupling_dir)
            
            print(f"✓ Successfully processed {coupling_name}")
            
        except Exception as e:
            print(f"✗ Error processing {coupling_name}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    if not successful_analyses:
        print("\n✗ No successful analyses to plot!")
        return
    
    print(f"\n{'='*70}")
    print(f"Successfully processed {len(successful_analyses)} coupling strengths")
    print(f"{'='*70}")
    
    # Generate plots
    try:
        print("\n" + "="*70)
        print("Generating unified material time collapse plot...")
        print("="*70)
        analyzer.plot_unified_material_time_collapse(successful_analyses, output_dir=str(base_dir))
        print("✓ Created: unified_material_time_collapse.png")
    except Exception as e:
        print(f"✗ Error creating unified collapse plot: {e}")
        import traceback
        traceback.print_exc()
    
    try:
        print("\n" + "="*70)
        print("Fitting stretched exponentials to collapsed data...")
        print("="*70)
        fit_results = analyzer.fit_stretched_exponential_to_collapse(successful_analyses, output_dir=str(base_dir))
        print("✓ Created: stretched_exponential_fits.png")
        print("✓ Created: stretched_exponential_fits.pdf")
    except Exception as e:
        print(f"✗ Error in stretched exponential fitting: {e}")
        import traceback
        traceback.print_exc()
    
    try:
        print("\n" + "="*70)
        print("Creating overlapped stretched exponential fits...")
        print("="*70)
        analyzer.plot_overlapped_stretched_exponential_fits(successful_analyses, output_dir=str(base_dir))
        print("✓ Created: overlapped_stretched_exponential_fits.png")
        print("✓ Created: overlapped_stretched_exponential_fits.pdf")
    except Exception as e:
        print(f"✗ Error in overlapped fits: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print("\nGenerated files:")
    print("  • unified_material_time_collapse.png")
    print("  • stretched_exponential_fits.png (+ .pdf)")
    print("  • overlapped_stretched_exponential_fits.png (+ .pdf)")
    print("  • stretched_exponential_fit_results.txt")
    print("\nAll plots use the CORRECTED material time reconstruction method!")
    print("="*80)


if __name__ == "__main__":
    import argparse

    from dataset_profile import add_profile_args
    from repro_bootstrap import setup_profile

    parser = argparse.ArgumentParser(
        description="Run corrected material time analysis and collapse plots"
    )
    add_profile_args(parser, default="paper")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path(
            "/home/mh7373/GitRepos/cav-hoomd/examples/final_nodiss_cavitymd"
        ),
        help="Output directory (paper archive layout) or override for profile runs",
    )
    args = parser.parse_args()
    profile = setup_profile(args, default="paper")
    if args.profile != "paper":
        output_dir = profile.staged_root
        if profile.figures_dir is not None:
            output_dir = profile.figures_dir
        run_corrected_analysis(output_dir, profile=profile)
    else:
        run_corrected_analysis(args.base_dir)

