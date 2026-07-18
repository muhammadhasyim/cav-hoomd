#!/usr/bin/env python3
"""
F(k,t) Analysis and Plotting Script

This script analyzes F(k,t) data from averaged replica files and creates three types of plots:
1. Multi-panel plot: F(k,t) vs time for different ref's (one panel per coupling strength)
2. Multi-panel plot: F(k,t) vs time for different coupling strengths (one panel per ref)
3. Two-panel plot: Relaxation times analysis

Usage: python plot_fkt_analysis.py [--base_dir ./]
"""

import argparse
import os
import sys
import glob
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for remote servers
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import re
from scipy.interpolate import interp1d
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import subprocess
from datetime import datetime

# Local KWW estimator (same directory when run as a script).
_SHARED_DIR = Path(__file__).resolve().parent
if str(_SHARED_DIR) not in sys.path:
    sys.path.insert(0, str(_SHARED_DIR))
from relaxation_fit import (  # noqa: E402
    DEFAULT_TARGET as KWW_TARGET,
    fit_kww_shared_beta,
    fit_kww_single,
    kww,
)

def parse_coupling_strength(folder_name):
    """Parse coupling strength from folder name."""
    coupling_part = folder_name.replace('cavity_coupling_', '')
    
    if '0epos00' in coupling_part:
        return 0.0, '0.0'
    elif 'eneg' in coupling_part:
        parts = coupling_part.split('eneg')
        if len(parts) == 2:
            mantissa = int(parts[0])
            exponent = int(parts[1][1])
            value = mantissa * (10 ** (-exponent))
            return value, f'{mantissa}\\times10^{{-{exponent}}}'
    
    return float('nan'), coupling_part

# When a dataset profile supplies the coupling catalogue directly (e.g. the
# in-repo aging_weak_lambda runs), the paper's fixed epsilon whitelist must be
# bypassed so all profile couplings are plotted.
_USE_PROFILE_COUPLINGS = False

# Axis symbol used in coupling-axis labels ('epsilon' for paper, 'lambda' for aging).
COUPLING_AXIS_LABEL = 'Coupling Strength'

# Parenthetical subtitle describing which coupling set is plotted.
COUPLING_SET_LABEL = 'Filtered: 0.0, 3×10⁻⁴, 5×10⁻⁴, 1×10⁻³'

# Math symbol for the coupling strength in legends/colorbars.
COUPLING_SYMBOL_TEX = 'g'
COUPLING_SYMBOL_PLAIN = 'g'


def filter_coupling_strengths(coupling_info):
    """Filter to only include specific coupling strengths: 0.0, 3e-4, 5e-4, 1e-3.

    When ``_USE_PROFILE_COUPLINGS`` is set (profile-driven runs), the whitelist is
    bypassed and every coupling with a finite value is retained.
    """
    if _USE_PROFILE_COUPLINGS:
        return {
            name: (value, label)
            for name, (value, label) in coupling_info.items()
            if not np.isnan(value)
        }

    target_couplings = [0.0, 3e-4, 5e-4, 7e-4, 1e-3, 2e-3, 4e-3, 1e-2]
    
    filtered_couplings = {}
    for coupling_name, (coupling_value, coupling_label) in coupling_info.items():
        if not np.isnan(coupling_value) and any(abs(coupling_value - target) < 1e-6 for target in target_couplings):
            filtered_couplings[coupling_name] = (coupling_value, coupling_label)
    
    return filtered_couplings


def collect_data_from_profile(profile):
    """Collect F(k,t) data using a dataset profile's staged coupling catalogue.

    Unlike :func:`collect_data`, coupling values and labels come directly from the
    profile (axis value + formatted label) instead of parsing epsilon tags from
    folder names, so aging (lambda-axis) datasets are handled correctly.

    Parameters
    ----------
    profile : DatasetProfile
        Loaded dataset profile providing staged coupling directories and axis values.

    Returns
    -------
    tuple(dict, dict)
        ``data_dict`` mapping folder name -> {ref_number: (time, fkt, counts)}
        and ``coupling_info`` mapping folder name -> (axis_value, label).
        ``counts`` may be None when the sample-count file is missing or misaligned.
    """
    data_dict = {}
    coupling_info = {}

    for entry in profile.couplings:
        folder = profile.staged_root / entry.staged_dir_name
        coupling_name = folder.name
        label = f'{entry.axis_value:g}'
        coupling_info[coupling_name] = (entry.axis_value, label)

        if not folder.is_dir():
            print(f"  Skipping {coupling_name}: staged directory missing")
            continue

        print(f"Processing folder: {coupling_name} ({profile.axis} = {label})")

        ref_files = list(folder.glob('master_fskt_ref*.txt'))
        if not ref_files:
            print(f"  No master_fskt_ref*.txt files found in {folder}")
            continue

        folder_data = {}
        for ref_file in ref_files:
            ref_match = re.search(r'master_fskt_ref_?(\d+)\.txt', ref_file.name)
            if not ref_match:
                continue
            ref_number = int(ref_match.group(1))
            time, fkt, counts = read_fkt_data(ref_file)
            if time is not None and fkt is not None:
                folder_data[ref_number] = (time, fkt, counts)

        if folder_data:
            data_dict[coupling_name] = folder_data
        print(f"  Total refs found: {len(folder_data)}")

    return data_dict, coupling_info

def process_fkt_data(time, fkt, normalization_value=None):
    """
    Process F(k,t) data: normalize and remove values below 0.001.
    Also auto-detect maximum time for plotting.
    """
    time = np.array(time, dtype=np.float64)
    fkt = np.array(fkt, dtype=np.float64)

    # Remove NaN values and zero values
    valid_mask = ~(np.isnan(time) | np.isnan(fkt)) & (fkt != 0)
    time_clean = time[valid_mask]
    fkt_clean = fkt[valid_mask]
    
    if len(time_clean) < 2:
        return None, None, None
    
    # Sort by time to ensure monotonic behavior
    sort_indices = np.argsort(time_clean)
    time_sorted = time_clean[sort_indices]
    fkt_sorted = fkt_clean[sort_indices]
    
    # Normalize using provided normalization value
    if normalization_value is not None and normalization_value != 0:
        fkt_normalized = fkt_sorted / normalization_value
    else:
        # Fallback: normalize by first value
        if len(fkt_sorted) > 0:
            fkt_normalized = fkt_sorted / fkt_sorted[0]
        else:
            return None, None, None
    
    # Remove values below 0.001 after normalization
    above_threshold_mask = fkt_normalized >= 0.001
    if not np.any(above_threshold_mask):
        return None, None, None
    
    time_filtered = time_sorted[above_threshold_mask]
    fkt_filtered = fkt_normalized[above_threshold_mask]
    
    # Auto-detect maximum time for plotting
    # Use the last time point where F(k,t) >= 0.001
    max_time = time_filtered[-1] if len(time_filtered) > 0 else None
    
    return time_filtered, fkt_filtered, max_time

def unpack_fkt_entry(entry):
    """Unpack a data_dict entry as ``(time, fkt, counts)``.

    Supports legacy 2-tuples ``(time, fkt)`` (counts=None) and 3-tuples
    ``(time, fkt, counts)``.
    """
    if entry is None:
        return None, None, None
    if len(entry) >= 3:
        return entry[0], entry[1], entry[2]
    return entry[0], entry[1], None


def read_sample_counts(counts_path, n_expected):
    """Read a one-integer-per-line sample-count file aligned to F(k,t) rows.

    Returns a float array of length ``n_expected``, or None if the file is
    missing or the length does not match (caller should then use uniform weights).
    """
    counts_path = Path(counts_path)
    if not counts_path.is_file():
        return None
    try:
        counts = np.loadtxt(counts_path, dtype=np.float64)
        if counts.ndim == 0:
            counts = np.array([float(counts)])
        if len(counts) != n_expected:
            print(
                f"  Warning: sample-count length mismatch for {counts_path.name}: "
                f"{len(counts)} vs {n_expected} data rows; using uniform weights"
            )
            return None
        return counts
    except Exception as e:
        print(f"  Warning: failed to read sample counts {counts_path}: {e}")
        return None


def sample_counts_path_for(ref_file):
    """Return the conventional sample-count path for a master F(k,t) file."""
    ref_file = Path(ref_file)
    # master_fskt_ref_001.txt -> fkt_sample_counts/master_fskt_ref_001_sample_counts.txt
    return ref_file.parent / "fkt_sample_counts" / f"{ref_file.stem}_sample_counts.txt"


def read_fkt_data(file_path, load_counts=True):
    """Read F(k,t) data from averaged replica file.

    Parameters
    ----------
    file_path : path-like
        Path to ``master_fskt_ref_XXX.txt``.
    load_counts : bool
        If True, also load the matching sample-count file when present and
        length-aligned.

    Returns
    -------
    tuple
        ``(time, fkt, counts)`` where ``counts`` may be None.
    """
    try:
        data = pd.read_csv(file_path, sep=r'\s+', comment='#')#, header=None)
        
        time = data.iloc[:, 0].values
        fkt = data.iloc[:, 1].values
        counts = None
        if load_counts:
            counts = read_sample_counts(sample_counts_path_for(file_path), len(time))
        return time, fkt, counts
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None, None, None


def _normalization_value_from_entry(entry):
    """First nonzero F(k,t) value from a ref0-style entry, or None."""
    _time, fkt, _counts = unpack_fkt_entry(entry)
    if fkt is None:
        return None
    nonzero_mask = fkt != 0
    if not np.any(nonzero_mask):
        return None
    return float(fkt[np.where(nonzero_mask)[0][0]])


def compute_paper_raw_relaxation_map(data_dict, coupling_names=None, target=KWW_TARGET):
    """Compute paper-style structural relaxation times for each (coupling, ref).

    For each reference curve, normalize by that curve's own first nonzero
    ``F`` value (proxy for :math:`S_k = F(k,0)`) and report the direct
    ``phi = target`` (default 0.1) crossing.  This matches arXiv:2603.15693
    Methods (``phi_k = F_k / S_k``, ``tau_s`` at ``phi_k = 0.1``) rather than
    a shared-beta KWW fit or a global ref0 normalization.

    Parameters
    ----------
    data_dict : dict
        ``coupling_name -> {ref_num: (time, fkt[, counts])}``.
    coupling_names : iterable or None
        Couplings to process; default is all keys in ``data_dict``.
    target : float
        Crossing level for ``tau_s``.

    Returns
    -------
    dict
        ``(coupling_name, ref_num) -> tau_s`` in picoseconds.
    """
    if coupling_names is None:
        coupling_names = list(data_dict.keys())

    results: dict[tuple[str, int], float] = {}
    for coupling_name in coupling_names:
        if coupling_name not in data_dict:
            continue
        folder_data = data_dict[coupling_name]
        for ref_num, entry in folder_data.items():
            time_arr, fkt_arr, _counts = unpack_fkt_entry(entry)
            if time_arr is None or fkt_arr is None:
                continue
            nonzero = fkt_arr != 0
            if not np.any(nonzero):
                continue
            norm_val = float(fkt_arr[np.where(nonzero)[0][0]])
            if norm_val == 0:
                continue
            tau = find_relaxation_time(
                time_arr, fkt_arr, target_value=target, normalization_value=norm_val
            )
            if np.isfinite(tau) and tau > 0:
                results[(coupling_name, int(ref_num))] = float(tau)
    return results


def compute_kww_relaxation_map(data_dict, coupling_names=None, target=KWW_TARGET):
    """Compute shared-beta KWW relaxation times for each (coupling, ref).

    For each coupling, fits all available reference curves jointly with one
    shared beta and per-ref ``(A, tau_K)``, then reports the analytic
    ``F = target`` (default 0.1) crossing of the fitted curve as ``tau_s``.

    Falls back to the raw single-crossing estimator when a KWW fit fails.

    Parameters
    ----------
    data_dict : dict
        ``coupling_name -> {ref_num: (time, fkt[, counts])}``.
    coupling_names : iterable or None
        Couplings to process; default is all keys in ``data_dict``.
    target : float
        Crossing level for ``tau_s``.

    Returns
    -------
    dict
        ``(coupling_name, ref_num) -> result_dict`` with keys
        ``tau_s``, ``A``, ``tau_K``, ``beta``, ``r_squared``, ``method``.
    """
    if coupling_names is None:
        coupling_names = list(data_dict.keys())

    results: dict[tuple[str, int], dict] = {}

    for coupling_name in coupling_names:
        if coupling_name not in data_dict:
            continue
        folder_data = data_dict[coupling_name]
        if not folder_data:
            continue

        norm_val = None
        if 0 in folder_data:
            norm_val = _normalization_value_from_entry(folder_data[0])

        ref_nums = sorted(folder_data.keys())
        curves = []
        meta = []  # (ref_num, time_raw, fkt_raw, counts)
        for ref_num in ref_nums:
            time_arr, fkt_arr, counts = unpack_fkt_entry(folder_data[ref_num])
            if time_arr is None or fkt_arr is None:
                continue
            nv = norm_val
            if nv is None or nv == 0:
                nonzero = fkt_arr != 0
                nv = float(fkt_arr[np.where(nonzero)[0][0]]) if np.any(nonzero) else None
            if nv is None or nv == 0:
                continue
            phi = np.asarray(fkt_arr, dtype=np.float64) / nv
            weights = counts if counts is not None else np.ones_like(phi)
            curves.append({
                "time": np.asarray(time_arr, dtype=np.float64),
                "phi": phi,
                "weights": np.asarray(weights, dtype=np.float64),
            })
            meta.append((ref_num, time_arr, fkt_arr, nv))

        if not curves:
            continue

        fit_list = fit_kww_shared_beta(curves, target=target)
        for j, (ref_num, time_arr, fkt_arr, nv) in enumerate(meta):
            fit = fit_list[j] if fit_list is not None else None
            if fit is not None and np.isfinite(fit.get("tau_s", np.nan)):
                results[(coupling_name, ref_num)] = {
                    "tau_s": float(fit["tau_s"]),
                    "A": float(fit["A"]),
                    "tau_K": float(fit["tau_K"]),
                    "beta": float(fit["beta"]),
                    "r_squared": float(fit["r_squared"]),
                    "method": "kww_shared_beta",
                }
            else:
                # Fallback: raw crossing at the same target.
                raw = find_relaxation_time(
                    time_arr, fkt_arr, target_value=target, normalization_value=nv
                )
                results[(coupling_name, ref_num)] = {
                    "tau_s": float(raw) if not np.isnan(raw) else float("nan"),
                    "A": float("nan"),
                    "tau_K": float("nan"),
                    "beta": float("nan"),
                    "r_squared": float("nan"),
                    "method": "raw_crossing_fallback",
                }

    return results

def find_relaxation_time(time, fkt, target_value=1.0 / np.e, normalization_value=None):
    """Find the relaxation time where F(k,t) = target_value.

    The default ``target_value`` is the 1/e criterion (F(k,t) = 1/e ~ 0.368),
    matching the labels used throughout this module.
    """
    try:
        # Remove NaN values and zero values (zero values are not meaningful for relaxation analysis)
        valid_mask = ~(np.isnan(time) | np.isnan(fkt)) & (fkt != 0)
        time_clean = time[valid_mask]
        fkt_clean = fkt[valid_mask]
        
        if len(time_clean) < 2:
            return np.nan
        
        # Sort by time to ensure monotonic behavior
        sort_indices = np.argsort(time_clean)
        time_sorted = time_clean[sort_indices]
        fkt_sorted = fkt_clean[sort_indices]
        
        # Don't artificially truncate - use all available non-zero data
        # Only remove last few points if they seem like noise (less aggressive)
        if len(fkt_sorted) > 10:
            n_keep = int(0.995 * len(fkt_sorted))  # Only remove 0.5% instead of 1%
            fkt_sorted = fkt_sorted[:n_keep]
            time_sorted = time_sorted[:n_keep]
        
        # Use provided normalization value
        if normalization_value is not None and normalization_value != 0:
            fkt_sorted = fkt_sorted / normalization_value
        else:
            # Fallback: normalize by first value (should not happen now since we have ref0 norm)
            if len(fkt_sorted) > 0:
                fkt_sorted = fkt_sorted / fkt_sorted[0]
            else:
                return np.nan
        
        # Check if target value is within the range of fkt
        if target_value < fkt_sorted.min() or target_value > fkt_sorted.max():
            return np.nan
        
        # Remove duplicate F(k,t) values by keeping the first occurrence
        _, unique_indices = np.unique(fkt_sorted, return_index=True)
        unique_indices = np.sort(unique_indices)
        fkt_unique = fkt_sorted[unique_indices]
        time_unique = time_sorted[unique_indices]
        
        if len(fkt_unique) < 2:
            return np.nan
        
        # Check again if target is in range after removing duplicates
        if target_value < fkt_unique.min() or target_value > fkt_unique.max():
            return np.nan
        
        # For F(k,t) data, we typically expect monotonic decay
        # If F(k,t) is not monotonic, we need to handle this carefully
        if not np.all(np.diff(fkt_unique) <= 0):  # Not monotonically decreasing
            # Find the part where F(k,t) crosses the target value (first crossing)
            crossing_indices = np.where((fkt_sorted[:-1] >= target_value) & (fkt_sorted[1:] <= target_value))[0]
            if len(crossing_indices) > 0:
                # Use linear interpolation for the first crossing
                idx = crossing_indices[0]
                t1, t2 = time_sorted[idx], time_sorted[idx+1]
                f1, f2 = fkt_sorted[idx], fkt_sorted[idx+1]
                
                if f1 != f2:  # Avoid division by zero
                    relaxation_time = t1 + (target_value - f1) * (t2 - t1) / (f2 - f1)
                    return relaxation_time if relaxation_time >= 0 else np.nan
                else:
                    return t1 if t1 >= 0 else np.nan
            else:
                return np.nan
        
        # If monotonic, use interpolation
        if len(fkt_unique) > 3:
            interp_func = interp1d(fkt_unique, time_unique, kind='cubic', 
                                 bounds_error=False, fill_value=np.nan)
        else:
            interp_func = interp1d(fkt_unique, time_unique, kind='linear',
                                 bounds_error=False, fill_value=np.nan)
        
        relaxation_time = interp_func(target_value)
        return relaxation_time if not np.isnan(relaxation_time) and relaxation_time >= 0 else np.nan
        
    except Exception as e:
        print(f"Error finding relaxation time: {e}")
        return np.nan

def collect_data(base_dir):
    """Collect all F(k,t) data from the coupling strength folders."""
    base_path = Path(base_dir)
    coupling_folders = list(base_path.glob('cavity_coupling_*200.0ps'))
    
    if not coupling_folders:
        print(f"No coupling folders found in {base_dir}")
        return {}, {}
    
    data_dict = {}
    coupling_info = {}
    
    for folder in coupling_folders:
        coupling_name = folder.name
        coupling_value, coupling_label = parse_coupling_strength(coupling_name)
        coupling_info[coupling_name] = (coupling_value, coupling_label)
        
        print(f"Processing folder: {coupling_name} (coupling = {coupling_label})")
        
        ref_files = list(folder.glob('master_fskt_ref*.txt'))
        
        if not ref_files:
            print(f"  No master_fskt_ref*.txt files found in {folder}")
            continue
        
        folder_data = {}
        
        for ref_file in ref_files:
            ref_match = re.search(r'master_fskt_ref_?(\d+)\.txt', ref_file.name)
            if ref_match:
                ref_number = int(ref_match.group(1))
                
                time, fkt, counts = read_fkt_data(ref_file)
                
                if time is not None and fkt is not None:
                    folder_data[ref_number] = (time, fkt, counts)
                    n_msg = f", {len(counts)} sample counts" if counts is not None else ", no sample counts"
                    print(f"  Read ref{ref_number}: {len(time)} data points{n_msg}")
                else:
                    print(f"  Failed to read ref{ref_number}")
        
        if folder_data:
            data_dict[coupling_name] = folder_data
        
        print(f"  Total refs found: {len(folder_data)}")
    
    return data_dict, coupling_info

def plot_fkt_by_coupling(data_dict, coupling_info, output_dir='.'):
    """Create multi-panel plot: F(k,t) vs time for different ref's."""
    print("\nCreating F(k,t) by coupling strength plot...")
    
    # Filter to specific coupling strengths
    filtered_coupling_info = filter_coupling_strengths(coupling_info)
    
    sorted_couplings = sorted(filtered_coupling_info.keys(), 
                            key=lambda x: filtered_coupling_info[x][0])
    
    n_couplings = len(sorted_couplings)
    if n_couplings == 0:
        print("No target coupling data to plot")
        return
    
    if n_couplings <= 2:
        rows, cols = 1, n_couplings
    elif n_couplings <= 4:
        rows, cols = 2, 2
    else:
        rows = int(np.ceil(n_couplings / 3))
        cols = 3
    
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows))
    if n_couplings == 1:
        axes = [axes]
    elif rows == 1 or cols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()
    
    # Track maximum time across all data for consistent x-axis
    global_max_time = 0
    
    for i, coupling_name in enumerate(sorted_couplings):
        if i >= len(axes):
            break
            
        ax = axes[i]
        
        # Check if this coupling has data
        if coupling_name not in data_dict:
            coupling_value, coupling_label = filtered_coupling_info[coupling_name]
            ax.text(0.5, 0.5, 'No data files found', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Coupling: {coupling_label}')
            continue
            
        folder_data = data_dict[coupling_name]
        coupling_value, coupling_label = filtered_coupling_info[coupling_name]
        
        ref_numbers = sorted(folder_data.keys())
        
        if not ref_numbers:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Coupling: {coupling_label}')
            continue
        
        # Get normalization value from ref0
        normalization_value = None
        if 0 in folder_data:
            time_ref0, fkt_ref0, _ = unpack_fkt_entry(folder_data[0])
            # Find first non-zero value from ref0 for normalization
            nonzero_mask = fkt_ref0 != 0
            if np.any(nonzero_mask):
                first_nonzero_idx = np.where(nonzero_mask)[0][0]
                normalization_value = fkt_ref0[first_nonzero_idx]
        
        norm = Normalize(vmin=0, vmax=max(ref_numbers))
        cmap = plt.colormaps.get_cmap('coolwarm')
        
        panel_max_time = 0
        
        for ref_num in ref_numbers:
            time, fkt, _ = unpack_fkt_entry(folder_data[ref_num])
            color = cmap(norm(ref_num))
            
            # Process F(k,t) data with filtering and normalization
            time_processed, fkt_processed, max_time = process_fkt_data(time, fkt, normalization_value)
            
            if time_processed is not None and fkt_processed is not None:
                ax.plot(time_processed, fkt_processed, color=color, linewidth=2, alpha=0.8, 
                       label=f'ref{ref_num}')
                
                if max_time is not None:
                    panel_max_time = max(panel_max_time, max_time)
                    global_max_time = max(global_max_time, max_time)
        
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('F(k,t) (normalized, ≥0.001)')
        ax.set_title(f'Coupling: {coupling_label}')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, ncol=2)
        ax.set_ylim(bottom=0.001)  # Start y-axis at threshold
        
        # Set x-axis limit based on data for this panel
        if panel_max_time > 0:
            ax.set_xlim(0, panel_max_time * 1.05)  # Add 5% padding
    
    for i in range(n_couplings, len(axes)):
        fig.delaxes(axes[i])
    
    plt.suptitle(f'F(k,t) vs Time: Different References by Coupling Strength\n({COUPLING_SET_LABEL})', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    output_file = Path(output_dir) / 'fkt_by_coupling_filtered.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_fkt_by_ref(data_dict, coupling_info, output_dir='.'):
    """Create multi-panel plot: F(k,t) vs time for different coupling strengths."""
    print("\nCreating F(k,t) by reference plot...")
    
    # Filter to specific coupling strengths
    filtered_coupling_info = filter_coupling_strengths(coupling_info)
    
    all_refs = set()
    for coupling_name in filtered_coupling_info.keys():
        if coupling_name in data_dict:
            all_refs.update(data_dict[coupling_name].keys())
    
    all_refs = sorted(all_refs)
    
    if not all_refs:
        print("No reference data to plot")
        return
    
    n_refs = len(all_refs)
    
    if n_refs <= 2:
        rows, cols = 1, n_refs
    elif n_refs <= 4:
        rows, cols = 2, 2
    else:
        rows = int(np.ceil(n_refs / 3))
        cols = 3
    
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows))
    if n_refs == 1:
        axes = [axes]
    elif rows == 1 or cols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()
    
    sorted_couplings = sorted(filtered_coupling_info.keys(), 
                            key=lambda x: filtered_coupling_info[x][0])
    
    coupling_values = [filtered_coupling_info[c][0] for c in sorted_couplings]
    norm = Normalize(vmin=min(coupling_values), vmax=max(coupling_values))
    cmap = plt.colormaps.get_cmap('coolwarm')
    
    # Track maximum time across all data for consistent x-axis
    global_max_time = 0
    
    for i, ref_num in enumerate(all_refs):
        if i >= len(axes):
            break
            
        ax = axes[i]
        has_data = False
        panel_max_time = 0
        
        for coupling_name in sorted_couplings:
            if coupling_name in data_dict and ref_num in data_dict[coupling_name]:
                time, fkt, _ = unpack_fkt_entry(data_dict[coupling_name][ref_num])
                coupling_value, coupling_label = filtered_coupling_info[coupling_name]
                
                color = cmap(norm(coupling_value))
                
                # Get normalization value from ref0 of this coupling
                normalization_value = None
                if 0 in data_dict[coupling_name]:
                    time_ref0, fkt_ref0, _ = unpack_fkt_entry(data_dict[coupling_name][0])
                    # Find first non-zero value from ref0 for normalization
                    nonzero_mask = fkt_ref0 != 0
                    if np.any(nonzero_mask):
                        first_nonzero_idx = np.where(nonzero_mask)[0][0]
                        normalization_value = fkt_ref0[first_nonzero_idx]
                
                # Process F(k,t) data with filtering and normalization
                time_processed, fkt_processed, max_time = process_fkt_data(time, fkt, normalization_value)
                
                if time_processed is not None and fkt_processed is not None:
                    ax.plot(time_processed, fkt_processed, color=color, linewidth=2, alpha=0.8, 
                       label=f'{coupling_label}')
                    has_data = True
                    
                    if max_time is not None:
                        panel_max_time = max(panel_max_time, max_time)
                        global_max_time = max(global_max_time, max_time)
        
        if not has_data:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('F(k,t) (normalized, ≥0.001)')
        ax.set_title(f'ref{ref_num}')
        ax.grid(True, alpha=0.3)
        
        if has_data:
            ax.legend(fontsize=8)
            ax.set_ylim(bottom=0.001)  # Start y-axis at threshold
            
            # Set x-axis limit based on data for this panel
            if panel_max_time > 0:
                ax.set_xlim(0, panel_max_time * 1.05)  # Add 5% padding
    
    for i in range(n_refs, len(axes)):
        fig.delaxes(axes[i])
    
    plt.suptitle(f'F(k,t) vs Time: Different Coupling Strengths by Reference\n({COUPLING_SET_LABEL})', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    output_file = Path(output_dir) / 'fkt_by_ref_filtered.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_fkt_diagnostic(data_dict, coupling_info, output_dir='.'):
    """Create diagnostic plot with KWW fit overlay and F=0.1 crossing marker."""
    print("\nCreating F(k,t) diagnostic plot...")
    
    # Filter to specific coupling strengths
    filtered_coupling_info = filter_coupling_strengths(coupling_info)
    
    sorted_couplings = sorted(filtered_coupling_info.keys(), 
                            key=lambda x: filtered_coupling_info[x][0])
    
    n_couplings = len(sorted_couplings)
    if n_couplings == 0:
        return

    kww_map = compute_kww_relaxation_map(
        data_dict, coupling_names=sorted_couplings[:4], target=KWW_TARGET
    )
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    target_line = KWW_TARGET  # F=0.1 crossing of the KWW fit
    global_max_time = 0
    
    for i, coupling_name in enumerate(sorted_couplings[:4]):  # Limit to 4 for 2x2 grid
        if i >= len(axes):
            break
            
        ax = axes[i]
        
        # Check if this coupling has data
        if coupling_name not in data_dict:
            coupling_value, coupling_label = filtered_coupling_info[coupling_name]
            ax.text(0.5, 0.5, 'No data files found', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Coupling: {coupling_label}')
            continue
            
        folder_data = data_dict[coupling_name]
        coupling_value, coupling_label = filtered_coupling_info[coupling_name]
        
        if 0 in folder_data:  # Only plot ref0 for diagnostic
            time, fkt, counts = unpack_fkt_entry(folder_data[0])
            fit = kww_map.get((coupling_name, 0))
            
            # Get normalization value from ref0 (same as current data since this is ref0)
            normalization_value = None
            nonzero_mask = fkt != 0
            if np.any(nonzero_mask):
                first_nonzero_idx = np.where(nonzero_mask)[0][0]
                normalization_value = fkt[first_nonzero_idx]
            
            # Process F(k,t) data with filtering and normalization
            time_processed, fkt_processed, max_time = process_fkt_data(time, fkt, normalization_value)
            
            if time_processed is not None and fkt_processed is not None:
                # Plot F(k,t) vs time
                ax.plot(time_processed, fkt_processed, 'b-', linewidth=2, alpha=0.8, label='F(k,t)')

                if fit is not None and np.isfinite(fit.get("tau_K", np.nan)):
                    t_fit = np.linspace(max(time_processed.min(), 1e-3),
                                        time_processed.max(), 400)
                    phi_fit = kww(t_fit, fit["A"], fit["tau_K"], fit["beta"])
                    ax.plot(t_fit, phi_fit, 'k--', linewidth=1.8, alpha=0.85,
                            label=f'KWW (β={fit["beta"]:.2f})')
                
                ax.axhline(y=target_line, color='red', linestyle='--', alpha=0.7, 
                      label=f'F = {target_line:.1f}')
                
                rel_time = fit["tau_s"] if fit is not None else np.nan
                if not np.isnan(rel_time) and rel_time > 0:
                    ax.axvline(x=rel_time, color='green', linestyle=':', alpha=0.7,
                              label=f'τ_s = {rel_time:.1f} ps')
                    ax.plot(rel_time, target_line, 'go', markersize=8, alpha=0.8)
                
                if max_time is not None:
                    global_max_time = max(global_max_time, max_time)
                    ax.set_xlim(0, max_time * 1.05)  # Add 5% padding

                fkt_min, fkt_max = fkt_processed.min(), fkt_processed.max()
                info_text = f'Range: {fkt_min:.3f} - {fkt_max:.3f}'
                if fit is not None and np.isfinite(fit.get("r_squared", np.nan)):
                    info_text += f'\nKWW R²={fit["r_squared"]:.3f}'
                ax.text(0.02, 0.98, info_text, transform=ax.transAxes, 
                       verticalalignment='top', fontsize=8,
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            else:
                # If all zeros, show message
                ax.text(0.5, 0.5, 'No valid data\n(all zeros or below threshold)', 
                       ha='center', va='center', transform=ax.transAxes)
            
            ax.set_xlabel('Time (ps)')
            ax.set_ylabel('F(k,t) (normalized, ≥0.001)')
            ax.set_title(f'Coupling: {coupling_label}')
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=9)
            ax.set_ylim(bottom=0.001)  # Start y-axis at threshold
        else:
            ax.text(0.5, 0.5, 'No ref0 data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Coupling: {coupling_label}')
    
    # Remove unused subplots
    for i in range(len(sorted_couplings), len(axes)):
        fig.delaxes(axes[i])
    
    plt.suptitle(
        f'F(k,t) Diagnostic: KWW Fit + F={KWW_TARGET} Crossing\n({COUPLING_SET_LABEL})',
        fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    output_file = Path(output_dir) / 'fkt_diagnostic_filtered.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def save_relaxation_analysis_data(data_dict, filtered_coupling_info, output_dir='.', kww_map=None):
    """Save relaxation analysis data to text files.

    Uses the shared-beta KWW estimator (F=0.1 crossing of the fitted curve).
    """
    print("\nSaving relaxation analysis data...")
    
    # Get sorted couplings and valid references
    sorted_couplings = sorted(filtered_coupling_info.keys(), 
                            key=lambda x: filtered_coupling_info[x][0])
    valid_couplings = [c for c in sorted_couplings if c in data_dict]
    
    # Determine all available references across all couplings
    all_refs = set()
    for coupling_name in valid_couplings:
        all_refs.update(data_dict[coupling_name].keys())
    all_refs = sorted(all_refs)

    switch_ps = 200.0

    if kww_map is None:
        kww_map = compute_kww_relaxation_map(
            data_dict, coupling_names=valid_couplings, target=KWW_TARGET
        )

    # Zero-coupling averaged baseline tau (matches the styled measured plot).
    ref_coupling_name = None
    for coupling_name in valid_couplings:
        if abs(filtered_coupling_info[coupling_name][0]) < 1e-10:
            ref_coupling_name = coupling_name
            break

    baseline_tau = None
    if ref_coupling_name is not None:
        baseline_vals = [
            kww_map[(ref_coupling_name, ref)]["tau_s"]
            for ref in all_refs
            if (ref_coupling_name, ref) in kww_map
            and np.isfinite(kww_map[(ref_coupling_name, ref)]["tau_s"])
            and kww_map[(ref_coupling_name, ref)]["tau_s"] > 0
        ]
        if baseline_vals:
            baseline_tau = float(np.mean(baseline_vals))

    if baseline_tau is None:
        print("  Warning: no valid zero-coupling baseline; tilde_tau_s will be NaN in exports")

    def _tilde(rel_time):
        if baseline_tau and not np.isnan(rel_time):
            return rel_time / baseline_tau
        return np.nan

    def _row_from_kww(coupling_name, ref_num):
        key = (coupling_name, ref_num)
        if key not in kww_map:
            return None
        fit = kww_map[key]
        rel_time = fit["tau_s"]
        if not np.isfinite(rel_time):
            return None
        normalization_value = None
        if 0 in data_dict[coupling_name]:
            normalization_value = _normalization_value_from_entry(data_dict[coupling_name][0])
        coupling_value, coupling_label = filtered_coupling_info[coupling_name]
        ref_time_ps = ref_num * switch_ps
        return {
            'ref_num': ref_num,
            'ref_time_ps': ref_time_ps,
            't_w_ps': ref_time_ps - switch_ps,
            'coupling_name': coupling_name,
            'coupling_value': coupling_value,
            'coupling_label': coupling_label,
            'relaxation_time_ps': rel_time,
            'tilde_tau_s': _tilde(rel_time),
            'A': fit.get("A", float("nan")),
            'tau_K': fit.get("tau_K", float("nan")),
            'beta': fit.get("beta", float("nan")),
            'r_squared': fit.get("r_squared", float("nan")),
            'method': fit.get("method", ""),
            'normalization_value': normalization_value if normalization_value is not None else float("nan"),
        }

    # Panel 1 Data: Relaxation time vs coupling strength for different references
    panel1_data = []
    for ref_num in all_refs:
        for coupling_name in valid_couplings:
            if ref_num in data_dict[coupling_name]:
                row = _row_from_kww(coupling_name, ref_num)
                if row is not None:
                    panel1_data.append(row)
    
    # Panel 2 Data: Relaxation time vs reference time for different coupling strengths
    panel2_data = []
    for coupling_name in valid_couplings:
        for ref_num in sorted(data_dict[coupling_name].keys()):
            row = _row_from_kww(coupling_name, ref_num)
            if row is not None:
                panel2_data.append(row)
    
    col_header = (
        f"{'ref_num':<8} {'ref_time_ps':<12} {'t_w_ps':<10} {'coupling_name':<40} "
        f"{'coupling_value':<15} {'coupling_label':<20} {'relaxation_time_ps':<18} "
        f"{'tilde_tau_s':<14} {'A':<10} {'tau_K':<14} {'beta':<10} {'r_squared':<10} "
        f"{'method':<22} {'normalization_value':<18}"
    )

    def _write_row(f, entry):
        f.write(
            f"{entry['ref_num']:<8} {entry['ref_time_ps']:<12.1f} {entry['t_w_ps']:<10.1f} "
            f"{entry['coupling_name']:<40} {entry['coupling_value']:<15.6e} "
            f"{entry['coupling_label']:<20} {entry['relaxation_time_ps']:<18.6f} "
            f"{entry['tilde_tau_s']:<14.6f} {entry['A']:<10.6f} {entry['tau_K']:<14.6f} "
            f"{entry['beta']:<10.6f} {entry['r_squared']:<10.6f} {entry['method']:<22} "
            f"{entry['normalization_value']:<18.6e}\n"
        )

    # Save Panel 1 data (Relaxation time vs coupling strength)
    panel1_file = Path(output_dir) / 'relaxation_vs_coupling_data.txt'
    with open(panel1_file, 'w') as f:
        f.write("# Relaxation Time vs Coupling Strength Data\n")
        f.write("# Generated from F(k,t) analysis: F=0.1 crossing of shared-beta KWW fit\n")
        f.write(f"# tilde_tau_s = relaxation_time_ps / baseline_tau; "
                f"baseline_tau = {baseline_tau if baseline_tau is not None else float('nan'):.6f} ps "
                f"(waiting-time-averaged zero-coupling tau)\n")
        f.write("# t_w_ps = ref_time_ps - 200 (waiting time referenced to cavity turn-on)\n")
        f.write("# Columns: ref_num, ref_time_ps, t_w_ps, coupling_name, coupling_value, "
                "coupling_label, relaxation_time_ps, tilde_tau_s, A, tau_K, beta, r_squared, "
                "method, normalization_value\n")
        f.write("#\n")
        f.write(col_header + "\n")
        for entry in panel1_data:
            _write_row(f, entry)
    
    print(f"  Saved Panel 1 data: {panel1_file}")
    
    # Save Panel 2 data (Relaxation time vs reference time)
    panel2_file = Path(output_dir) / 'relaxation_vs_reference_time_data.txt'
    # Same columns as panel1 for consistency (order differs only in row ordering).
    with open(panel2_file, 'w') as f:
        f.write("# Relaxation Time vs Reference Time Data\n")
        f.write("# Generated from F(k,t) analysis: F=0.1 crossing of shared-beta KWW fit\n")
        f.write(f"# tilde_tau_s = relaxation_time_ps / baseline_tau; "
                f"baseline_tau = {baseline_tau if baseline_tau is not None else float('nan'):.6f} ps "
                f"(waiting-time-averaged zero-coupling tau)\n")
        f.write("# t_w_ps = ref_time_ps - 200 (waiting time referenced to cavity turn-on)\n")
        f.write("# Columns: ref_num, ref_time_ps, t_w_ps, coupling_name, coupling_value, "
                "coupling_label, relaxation_time_ps, tilde_tau_s, A, tau_K, beta, r_squared, "
                "method, normalization_value\n")
        f.write("#\n")
        f.write(col_header + "\n")
        for entry in panel2_data:
            _write_row(f, entry)
    
    print(f"  Saved Panel 2 data: {panel2_file}")
    
    # Save combined summary data
    summary_file = Path(output_dir) / 'relaxation_analysis_summary.txt'
    with open(summary_file, 'w') as f:
        f.write("# Relaxation Analysis Summary\n")
        f.write("# Generated from F(k,t) analysis: F=0.1 crossing of shared-beta KWW fit\n")
        f.write(f"# Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# tilde_tau_s baseline (waiting-time-averaged zero-coupling tau): "
                f"{baseline_tau if baseline_tau is not None else float('nan'):.6f} ps\n")
        f.write("#\n")
        f.write("# PANEL 1: Relaxation Time vs Coupling Strength\n")
        f.write(f"# Total data points: {len(panel1_data)}\n")
        f.write(f"# Reference times analyzed: {sorted(set(entry['ref_time_ps'] for entry in panel1_data))}\n")
        f.write(f"# Coupling strengths: {sorted(set(entry['coupling_value'] for entry in panel1_data))}\n")
        f.write("#\n")
        f.write("# PANEL 2: Relaxation Time vs Reference Time\n")
        f.write(f"# Total data points: {len(panel2_data)}\n")
        f.write(f"# Coupling strengths analyzed: {sorted(set(entry['coupling_value'] for entry in panel2_data))}\n")
        if panel2_data:
            f.write(f"# Reference time range: {min(entry['ref_time_ps'] for entry in panel2_data):.1f} - {max(entry['ref_time_ps'] for entry in panel2_data):.1f} ps\n")
        f.write("#\n")
        f.write("# Files generated:\n")
        f.write(f"# - {panel1_file.name}: Panel 1 data (relaxation time vs coupling strength)\n")
        f.write(f"# - {panel2_file.name}: Panel 2 data (relaxation time vs reference time)\n")
        f.write(f"# - relaxation_analysis_filtered.png: Plot visualization\n")
    
    print(f"  Saved summary: {summary_file}")


def _latex_usable():
    """Return True only if a LaTeX toolchain can actually compile a trivial doc.

    A bare ``latex --version`` succeeds even when the format files are missing
    (broken TeX install), which later crashes matplotlib's usetex renderer. This
    performs a real minimal compile so broken installs fall back to mathtext.
    """
    import tempfile
    latex = shutil.which('latex')
    if latex is None:
        return False
    try:
        with tempfile.TemporaryDirectory() as tmp:
            tex = os.path.join(tmp, 'probe.tex')
            with open(tex, 'w', encoding='utf-8') as fh:
                fh.write(r'\documentclass{article}\begin{document}$x$\end{document}')
            subprocess.run(
                [latex, '-interaction=nonstopmode', '-halt-on-error', 'probe.tex'],
                cwd=tmp, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                timeout=30, check=True,
            )
        return True
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError):
        return False


def plot_measured_structural_relaxation_styled(
    data_dict,
    coupling_info,
    output_dir='.',
    fkt_kmag=None,
    include_kww_diagnostic=False,
):
    """Create paper-style two-panel plot of measured structural relaxation time.

    Left panel: normalized :math:`\\tilde{\\tau}_\\mathrm{s}` vs coupling strength,
    colored by waiting time :math:`t_\\mathrm{w}` (referenced to cavity turn-on).

    Right panel: :math:`\\tilde{\\tau}_\\mathrm{s}` vs :math:`t_\\mathrm{w}` for
    each coupling strength.  Data come from replica-averaged :math:`F(k,t)` masters.

    Relaxation times use the paper definition: per-curve
    :math:`\\phi = F / F(t\\to 0)` and the direct :math:`\\phi = 0.1` crossing
    (arXiv:2603.15693).  A shared-beta KWW estimator is optional diagnostics only.

    Normalization: each :math:`\\tau` is divided by the waiting-time-averaged
    zero-coupling relaxation time.  Because the zero-coupling run has no cavity and
    stays at equilibrium, its :math:`\\tau` is (statistically) independent of
    :math:`t_\\mathrm{w}`; averaging over all references collapses it to a single
    constant baseline instead of a per-reference value, which avoids injecting the
    baseline's uncorrelated statistical noise into :math:`\\tilde{\\tau}_\\mathrm{s}`.
    """
    print("\nCreating measured structural relaxation analysis (paper-style panels)...")

    filtered_coupling_info = filter_coupling_strengths(coupling_info)
    switch_ps = 200.0

    ref_coupling_name = None
    for name, (val, _label) in filtered_coupling_info.items():
        if abs(val) < 1e-10:
            ref_coupling_name = name
            break
    if ref_coupling_name is None or ref_coupling_name not in data_dict:
        print("  Warning: missing zero-coupling reference; skipping styled measured plot")
        return

    plt.style.use('classic')
    use_latex = _latex_usable()
    if use_latex:
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
    else:
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family'] = 'DejaVu Sans'

    def fmt(tex, fallback):
        return tex if use_latex else fallback

    raw_taus = compute_paper_raw_relaxation_map(
        data_dict, coupling_names=list(filtered_coupling_info.keys()), target=KWW_TARGET
    )
    print(f"  Paper-style direct phi=0.1 crossings: {len(raw_taus)}")
    if include_kww_diagnostic:
        kww_map = compute_kww_relaxation_map(
            data_dict, coupling_names=list(filtered_coupling_info.keys()), target=KWW_TARGET
        )
        n_kww = sum(1 for v in kww_map.values() if v.get("method") == "kww_shared_beta")
        print(f"  KWW diagnostic fits (not plotted): {n_kww}")

    all_refs = sorted({ref for _name, ref in raw_taus})
    waiting_times_lab = [ref * switch_ps for ref in all_refs if ref * switch_ps > 0]

    # Zero-coupling baseline: average tau over all waiting times. The reference run
    # has no cavity and remains at equilibrium, so its tau should be independent of
    # t_w; collapsing it to a single constant avoids injecting the baseline's
    # uncorrelated statistical scatter into the normalized tilde{tau}_s.
    ref_baseline_taus = [
        raw_taus[(ref_coupling_name, ref)]
        for ref in all_refs
        if (ref_coupling_name, ref) in raw_taus and raw_taus[(ref_coupling_name, ref)] > 0
    ]
    if not ref_baseline_taus:
        print("  Warning: no valid zero-coupling relaxation times; skipping styled measured plot")
        return
    ref_tau_baseline = float(np.mean(ref_baseline_taus))
    print(f"  Zero-coupling baseline tau_s (averaged): {ref_tau_baseline:.3f} ps")
    if fkt_kmag is not None:
        print(
            f"  Note: F(k,t) masters were collected at |k|={fkt_kmag:g} "
            f"(paper Fig. 2 uses |k|=6.02)."
        )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    waiting_time_norm = Normalize(vmin=0, vmax=max((tw - switch_ps for tw in waiting_times_lab), default=1))
    waiting_time_cmap = plt.colormaps.get_cmap('viridis')

    for tw_lab in waiting_times_lab:
        ref_num = int(round(tw_lab / switch_ps))
        tw_shifted = tw_lab - switch_ps
        xs, ys = [], []
        for coupling_name, (coupling_value, _label) in sorted(
            filtered_coupling_info.items(), key=lambda item: item[1][0]
        ):
            key = (coupling_name, ref_num)
            if key not in raw_taus:
                continue
            xs.append(coupling_value)
            ys.append(raw_taus[key] / ref_tau_baseline)
        if not xs:
            continue
        color = waiting_time_cmap(waiting_time_norm(tw_shifted))
        ax1.plot(xs, ys, '-D', color=color, linewidth=3, markersize=8,
                 markerfacecolor='white', markeredgecolor=color, markeredgewidth=2.0)

    ax1.set_xlabel(fmt(r'$\lambda$ (a.u.)' if COUPLING_SYMBOL_TEX == r'\lambda'
                       else rf'${COUPLING_SYMBOL_TEX}$ (a.u.)',
                       'λ (a.u.)' if COUPLING_SYMBOL_PLAIN == 'λ' else f'{COUPLING_SYMBOL_PLAIN} (a.u.)'),
               fontsize=16)
    ax1.set_ylabel(fmt(r'$\tilde{\tau}_{\mathrm{s}}$', 'τ̃_s'), fontsize=16)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.axhline(y=1, color='gray', linestyle='--', alpha=0.7, linewidth=1.5)
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0.7)
    sm1 = cm.ScalarMappable(norm=waiting_time_norm, cmap=waiting_time_cmap)
    sm1.set_array([])
    cbar1 = plt.colorbar(sm1, ax=ax1, location='top', shrink=0.8, pad=0.1)
    cbar1.set_label(fmt(r'$t_{\mathrm{w}}$ (ps)', 't_w (ps)'), fontsize=12, labelpad=10)

    coupling_values_sorted = sorted(
        filtered_coupling_info[name][0] for name in filtered_coupling_info if name in data_dict
    )
    coupling_norm = Normalize(vmin=0, vmax=max(coupling_values_sorted, default=1))
    coupling_cmap = plt.colormaps.get_cmap('coolwarm')

    for coupling_name, (coupling_value, coupling_label) in sorted(
        filtered_coupling_info.items(), key=lambda item: item[1][0]
    ):
        if coupling_name not in data_dict:
            continue
        xs, ys = [], []
        for ref_num in sorted(data_dict[coupling_name]):
            tw_lab = ref_num * switch_ps
            if tw_lab <= 0:
                continue
            key = (coupling_name, ref_num)
            if key not in raw_taus:
                continue
            xs.append(tw_lab - switch_ps)
            ys.append(raw_taus[key] / ref_tau_baseline)
        if not xs:
            continue
        color = coupling_cmap(coupling_norm(coupling_value))
        label = fmt(rf'$\lambda = {coupling_label}$', f'λ = {coupling_label}')
        ax2.plot(xs, ys, '-D', color=color, linewidth=3, markersize=8,
                 markerfacecolor='white', markeredgecolor=color, markeredgewidth=2.0,
                 label=label)

    ax2.set_xlabel(fmt(r'$t_{\mathrm{w}}$ (ps)', 't_w (ps)'), fontsize=16)
    ax2.set_ylabel(fmt(r'$\tilde{\tau}_{\mathrm{s}}$', 'τ̃_s'), fontsize=16)
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.axhline(y=1, color='gray', linestyle='--', alpha=0.7, linewidth=1.5)
    ax2.legend(fontsize=10, loc='best')
    ax2.set_xlim(left=0)
    ax2.set_ylim(bottom=0.7)

    if fkt_kmag is not None:
        fig.suptitle(
            fmt(
                rf'Direct $\phi=0.1$ crossing; data $|k|={fkt_kmag:g}$ '
                r'(paper $|k|=6.02$)',
                f'Direct φ=0.1 crossing; data |k|={fkt_kmag:g} (paper |k|=6.02)',
            ),
            fontsize=11,
            y=1.02,
        )

    plt.tight_layout(pad=2.0)
    out_base = Path(output_dir) / 'measured_structural_relaxation_analysis'
    plt.savefig(out_base.with_suffix('.png'), dpi=300, bbox_inches='tight')
    plt.savefig(out_base.with_suffix('.pdf'), dpi=300, bbox_inches='tight')
    print(f"Saved: {out_base.with_suffix('.png')}")
    print(f"Saved: {out_base.with_suffix('.pdf')}")
    plt.close()
    plt.rcParams.update(plt.rcParamsDefault)


def plot_relaxation_analysis(data_dict, coupling_info, output_dir='.'):
    """Create two-panel plot: Relaxation time analysis."""
    print("\nCreating relaxation time analysis plot...")
    
    # Set matplotlib style to classic and try to enable LaTeX
    plt.style.use('classic')
    if _latex_usable():
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        use_latex = True
        print("  Using LaTeX for mathematical typesetting")
    else:
        print("  Warning: LaTeX not available/usable, using default matplotlib rendering")
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family'] = 'DejaVu Sans'
        use_latex = False
    
    # Helper function for LaTeX formatting
    def latex_format(text, fallback_text=None):
        if use_latex:
            return text
        else:
            return fallback_text if fallback_text else text.replace('$', '').replace(r'\tau', 'τ').replace(r'\times', '×')
    
    # Filter to specific coupling strengths
    filtered_coupling_info = filter_coupling_strengths(coupling_info)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    
    sorted_couplings = sorted(filtered_coupling_info.keys(), 
                            key=lambda x: filtered_coupling_info[x][0])
    
    # Determine all available references across all couplings
    all_refs = set()
    for coupling_name in sorted_couplings:
        if coupling_name in data_dict:
            all_refs.update(data_dict[coupling_name].keys())
    all_refs = sorted(all_refs)
    
    # Panel 1: Relaxation time vs coupling strength for different references
    print("  Calculating relaxation times vs coupling strength for different references...")
    
    coupling_values = [filtered_coupling_info[c][0] for c in sorted_couplings if c in data_dict]
    coupling_labels = [filtered_coupling_info[c][1] for c in sorted_couplings if c in data_dict]
    
    valid_couplings = [c for c in sorted_couplings if c in data_dict]

    kww_map = compute_kww_relaxation_map(
        data_dict, coupling_names=valid_couplings, target=KWW_TARGET
    )
    
    # Create color mapping for reference times (Panel 1)
    ref_times_ps = [ref_num * 200 for ref_num in all_refs]
    ref_time_norm = Normalize(vmin=0, vmax=max(ref_times_ps) if ref_times_ps else 1)
    ref_time_cmap = plt.colormaps.get_cmap('viridis')
    
    # Store data for Panel 1 colorbar
    panel1_scatter_points = []
    
    for i, ref_num in enumerate(all_refs):
        ref_relaxation_times = []
        ref_coupling_values = []
        ref_coupling_labels = []
        
        for coupling_name in valid_couplings:
            key = (coupling_name, ref_num)
            if key in kww_map and np.isfinite(kww_map[key]["tau_s"]):
                coupling_value, coupling_label = filtered_coupling_info[coupling_name]
                ref_relaxation_times.append(kww_map[key]["tau_s"])
                ref_coupling_values.append(coupling_value)
                ref_coupling_labels.append(coupling_label)
        
        if ref_relaxation_times:  # Only plot if we have data
            ref_time_ps = ref_num * 200
            color = ref_time_cmap(ref_time_norm(ref_time_ps))
            
            scatter = ax1.scatter(ref_coupling_values, ref_relaxation_times, 
                               c=[ref_time_ps]*len(ref_relaxation_times), 
                               cmap='viridis', norm=ref_time_norm,
                               s=80, alpha=0.8, edgecolors='black', linewidth=0.5,
                               label=latex_format(f'$t = {ref_time_ps}$ ps', f't = {ref_time_ps} ps'))
            
            # Connect points with lines
            ax1.plot(ref_coupling_values, ref_relaxation_times, '-', 
                    color=color, linewidth=2, alpha=0.6)
            
            panel1_scatter_points.append(scatter)
            
            print(f"    ref{ref_num} (t={ref_time_ps} ps): {len(ref_relaxation_times)} valid coupling points")
    
    ax1.set_xlabel(latex_format(COUPLING_AXIS_LABEL, COUPLING_AXIS_LABEL.replace('$', '').replace(r'\lambda', 'λ')))
    ax1.set_ylabel(latex_format(r'Relaxation Time $\tau$ (ps)', 'Relaxation Time τ (ps)'))
    ax1.set_title(latex_format(
        r'Relaxation Time vs Coupling Strength' + '\n' + r'($F=0.1$ crossing of KWW fit)',
        'Relaxation Time vs Coupling Strength\n(F=0.1 crossing of KWW fit)'))
    ax1.grid(True, alpha=0.3)
    
    if coupling_values:
        ax1.set_xticks(coupling_values)
        if use_latex:
            ax1.set_xticklabels([r'$' + label + r'$' for label in coupling_labels], rotation=45)
        else:
            ax1.set_xticklabels(coupling_labels, rotation=45)
    
    # Set x-axis to start at 0 for Panel 1
    if coupling_values:
        ax1.set_xlim(left=0)
    
    # Add colorbar for Panel 1 (reference times)
    if panel1_scatter_points:
        cbar1 = plt.colorbar(panel1_scatter_points[0], ax=ax1, shrink=0.8)
        cbar1.set_label(latex_format(r'Reference Time $t$ (ps)', 'Reference Time t (ps)'), rotation=270, labelpad=20)
    
    # Panel 2: Relaxation time vs reference time for different coupling strengths
    print("  Calculating relaxation times vs reference time for different couplings...")
    
    # Create color mapping for coupling strengths (Panel 2)
    coupling_values_sorted = sorted([filtered_coupling_info[c][0] for c in valid_couplings])
    coupling_norm = Normalize(vmin=0, vmax=max(coupling_values_sorted))
    coupling_cmap = plt.colormaps.get_cmap('coolwarm')
    
    # Store data for Panel 2 colorbar
    panel2_scatter_points = []
    
    for i, coupling_name in enumerate(valid_couplings):
        coupling_value, coupling_label = filtered_coupling_info[coupling_name]
        folder_data = data_dict[coupling_name]
        
        ref_numbers = sorted(folder_data.keys())
        ref_times = []
        relaxation_times = []
        
        for ref_num in ref_numbers:
            key = (coupling_name, ref_num)
            if key in kww_map and np.isfinite(kww_map[key]["tau_s"]):
                ref_times.append(ref_num * 200)
                relaxation_times.append(kww_map[key]["tau_s"])
        
        if relaxation_times:  # Only plot if we have data
            color = coupling_cmap(coupling_norm(coupling_value))
            
            scatter = ax2.scatter(ref_times, relaxation_times,
                               c=[coupling_value]*len(relaxation_times),
                               cmap='coolwarm', norm=coupling_norm,
                               s=80, alpha=0.8, edgecolors='black', linewidth=0.5,
                               label=latex_format(f'${COUPLING_SYMBOL_TEX} = {coupling_label}$', 
                                                 f'{COUPLING_SYMBOL_PLAIN} = {coupling_label}'))
            
            # Connect points with lines
            ax2.plot(ref_times, relaxation_times, '-',
                    color=color, linewidth=2, alpha=0.6)
            
            panel2_scatter_points.append(scatter)
            
            print(f"    Coupling {coupling_label}: {len(relaxation_times)} valid reference points")
    
    # Highlight coupling turn-on time at 200 ps
    ax2.axvline(x=200, color='red', linestyle='--', linewidth=2, alpha=0.8,
                label=latex_format(r'Coupling Turn-On ($t = 200$ ps)', 'Coupling Turn-On (t = 200 ps)'))
    
    # Add shaded region to emphasize before/after coupling
    if ax2.get_xlim()[0] < 200:
        ax2.axvspan(ax2.get_xlim()[0], 200, alpha=0.1, color='gray', 
                   label=latex_format(r'Pre-Coupling Phase', 'Pre-Coupling Phase'))
    
    ax2.set_xlabel(latex_format(r'Reference Time $t$ (ps)', 'Reference Time t (ps)'))
    ax2.set_ylabel(latex_format(r'Relaxation Time $\tau$ (ps)', 'Relaxation Time τ (ps)'))
    ax2.set_title(latex_format(
        r'Relaxation Time vs Reference Time' + '\n' + r'($F=0.1$ crossing of KWW fit)',
        'Relaxation Time vs Reference Time\n(F=0.1 crossing of KWW fit)'))
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9, loc='best')
    
    # Set x-axis to start at 0 for Panel 2
    ax2.set_xlim(left=0)
    
    # Add colorbar for Panel 2 (coupling strengths)
    if panel2_scatter_points:
        cbar2 = plt.colorbar(panel2_scatter_points[0], ax=ax2, shrink=0.8)
        cbar2.set_label(latex_format(rf'Coupling Strength ${COUPLING_SYMBOL_TEX}$', f'Coupling Strength {COUPLING_SYMBOL_PLAIN}'), rotation=270, labelpad=20)
        
        # Format colorbar ticks for coupling strengths
        coupling_ticks = coupling_values_sorted
        if use_latex:
            coupling_tick_labels = [filtered_coupling_info[c][1].replace('×', r'$\times$') 
                                   for c in valid_couplings]
            cbar2.set_ticks(coupling_ticks)
            cbar2.set_ticklabels([f'${label}$' for label in coupling_tick_labels])
        else:
            coupling_tick_labels = [filtered_coupling_info[c][1] for c in valid_couplings]
            cbar2.set_ticks(coupling_ticks)
            cbar2.set_ticklabels(coupling_tick_labels)
    
    plt.suptitle(latex_format(
        r'Relaxation Time Analysis: $F=0.1$ Crossing of KWW Fit' + '\n' +
        r'Cavity MD with Coupling Turn-On at $t = 200$ ps',
        'Relaxation Time Analysis: F=0.1 Crossing of KWW Fit\nCavity MD with Coupling Turn-On at t = 200 ps'),
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    output_file = Path(output_dir) / 'relaxation_analysis_filtered.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()
    
    # Reset matplotlib settings
    plt.rcParams.update(plt.rcParamsDefault)
    
    # Save the data to a text file (reuse the same KWW map)
    save_relaxation_analysis_data(
        data_dict, filtered_coupling_info, output_dir, kww_map=kww_map
    )

def main():
    """Main function."""
    global _USE_PROFILE_COUPLINGS, COUPLING_AXIS_LABEL, COUPLING_SET_LABEL
    global COUPLING_SYMBOL_TEX, COUPLING_SYMBOL_PLAIN
    parser = argparse.ArgumentParser(description='Analyze and plot F(k,t) data from averaged replica files')
    parser.add_argument('--base_dir', default='./', 
                       help='Base directory containing coupling folders (default: ./)')
    parser.add_argument('--output_dir', default='./', 
                       help='Output directory for plots (default: ./)')
    parser.add_argument('--profile', default=None,
                       help='Dataset profile name; when set, couplings come from the profile '
                            'catalogue (bypasses the paper epsilon whitelist)')
    parser.add_argument('--staged-root', type=Path, default=None,
                       help='Override staged data root from the profile')
    parser.add_argument(
        '--fkt-kmag',
        type=float,
        default=None,
        help='Wavevector magnitude used when F(k,t) was collected '
             '(annotated on the measured tau plot; paper uses 6.02)',
    )
    
    args = parser.parse_args()
    
    # Set matplotlib style to classic globally
    plt.style.use('classic')
    
    print("F(k,t) Analysis and Plotting Script")
    print("=" * 50)

    if args.profile:
        REPRO_ROOT = Path(__file__).resolve().parents[1]
        sys.path.insert(0, str(REPRO_ROOT / "shared"))
        from repro_bootstrap import setup_profile

        profile = setup_profile(args, default=args.profile)
        _USE_PROFILE_COUPLINGS = True
        symbol = 'λ' if profile.axis == 'lambda' else 'ε'
        values = ', '.join(f'{e.axis_value:g}' for e in profile.couplings)
        COUPLING_SET_LABEL = f'{symbol} = {values}'
        if profile.axis == 'lambda':
            COUPLING_AXIS_LABEL = r'Coupling Strength $\lambda$'
            COUPLING_SYMBOL_TEX = r'\lambda'
            COUPLING_SYMBOL_PLAIN = 'λ'
        print(f"Profile: {profile.name} (axis={profile.axis})")
        print(f"Staged root: {profile.staged_root}")
        print(f"Output directory: {Path(args.output_dir).resolve()}")
        print("\nCollecting F(k,t) data...")
        data_dict, coupling_info = collect_data_from_profile(profile)
    else:
        print(f"Base directory: {Path(args.base_dir).resolve()}")
        print(f"Output directory: {Path(args.output_dir).resolve()}")
        print("\nCollecting F(k,t) data...")
        data_dict, coupling_info = collect_data(args.base_dir)
    
    if not data_dict:
        print("No data found! Make sure you have cavity_coupling_* folders with master_fskt_ref*.txt files.")
        return 1
    
    print(f"\nSummary:")
    print(f"  Found {len(data_dict)} coupling strengths")
    
    total_refs = set()
    for folder_data in data_dict.values():
        total_refs.update(folder_data.keys())
    
    print(f"  Total unique references: {sorted(total_refs)}")
    
    # Create output directory if it doesn't exist
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    # Generate all plots
    plot_fkt_by_coupling(data_dict, coupling_info, args.output_dir)
    plot_fkt_by_ref(data_dict, coupling_info, args.output_dir)
    plot_fkt_diagnostic(data_dict, coupling_info, args.output_dir)
    plot_relaxation_analysis(data_dict, coupling_info, args.output_dir)
    fkt_kmag = args.fkt_kmag
    if fkt_kmag is None and args.profile == "aging_weak_lambda":
        # Legacy campaign collected F(k,t) at |k|=1.0085 (paper uses 6.02).
        fkt_kmag = 1.0085
    plot_measured_structural_relaxation_styled(
        data_dict,
        coupling_info,
        args.output_dir,
        fkt_kmag=fkt_kmag,
    )
    
    print("\n" + "=" * 50)
    print("Analysis complete! Generated plots:")
    print(f"  1. fkt_by_coupling_filtered.png - F(k,t) for different refs (by coupling)")
    print(f"  2. fkt_by_ref_filtered.png - F(k,t) for different couplings (by ref)")
    print(f"  3. fkt_diagnostic_filtered.png - F(k,t) diagnostic with KWW fit + F=0.1 crossing")
    print(f"  4. relaxation_analysis_filtered.png - Relaxation time analysis (KWW F=0.1)")
    print(
        "  5. measured_structural_relaxation_analysis.png - "
        "Measured τ̃_s vs λ and t_w (paper direct φ=0.1)"
    )
    
    return 0

if __name__ == "__main__":
    exit(main())
