#!/usr/bin/env python3
"""
Finite-Size Scaling Study - Results Collection Script

This script collects results from all completed replicas and merges them
into a single summary CSV file. It also validates that all expected
replicas have completed successfully.

Usage:
    python collect_results.py --output-dir /path/to/output

Options:
    --output-dir DIR    Directory containing simulation outputs
    --summary-file FILE Output CSV file name (default: full_summary.csv)
    --validate          Check for missing/failed replicas
    --verbose           Print detailed progress
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# Expected replica counts for each system size
EXPECTED_REPLICAS = {
    100: 1250,
    250: 500,
    500: 250,
    1000: 125,
    2500: 50,
    5000: 25,
    10000: 13,
}


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Collect results from finite-size scaling study."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Directory containing simulation outputs",
    )
    parser.add_argument(
        "--summary-file",
        type=str,
        default="full_summary.csv",
        help="Output CSV file name (default: full_summary.csv)",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Check for missing/failed replicas",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed progress",
    )
    return parser.parse_args()


def get_replica_status(replica_dir: Path) -> str:
    """Read status file from replica directory."""
    status_file = replica_dir / "status.txt"
    if status_file.exists():
        return status_file.read_text().strip()
    return "unknown"


def find_replica_csv(replica_dir: Path) -> Path | None:
    """Find the summary CSV file in a replica directory."""
    # Look for CSV files matching expected patterns
    patterns = [
        "finite_size_scaling_summary.csv",
        "*_summary.csv",
        "*.csv",
    ]
    for pattern in patterns:
        matches = list(replica_dir.glob(pattern))
        if matches:
            return matches[0]
    return None


def collect_replica_data(
    output_dir: Path,
    n_molecules: int,
    replica_id: int,
    verbose: bool = False,
) -> Dict | None:
    """Collect data from a single replica."""
    replica_dir = output_dir / f"N_{n_molecules}" / f"replica_{replica_id}"
    
    if not replica_dir.exists():
        if verbose:
            print(f"  Missing: N={n_molecules}, replica={replica_id}")
        return None
    
    status = get_replica_status(replica_dir)
    
    if status != "success":
        if verbose:
            print(f"  Failed ({status}): N={n_molecules}, replica={replica_id}")
        return {
            "N_molecules": n_molecules,
            "replica": replica_id,
            "status": status,
            "lambda_single": None,
            "collective_coupling": None,
            "runtime_ps": None,
        }
    
    # Try to find and parse the CSV file
    csv_file = find_replica_csv(replica_dir)
    if csv_file and csv_file.exists():
        try:
            with csv_file.open("r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    # Return the first (and should be only) data row
                    return {
                        "N_molecules": int(row.get("N_molecules", n_molecules)),
                        "replica": int(row.get("replica", replica_id)),
                        "status": row.get("status", status),
                        "lambda_single": row.get("lambda_single"),
                        "collective_coupling": row.get("collective_coupling"),
                        "runtime_ps": row.get("runtime_ps"),
                        "device": row.get("device"),
                    }
        except Exception as e:
            if verbose:
                print(f"  Error reading CSV for N={n_molecules}, replica={replica_id}: {e}")
    
    # If no CSV found, return basic info
    return {
        "N_molecules": n_molecules,
        "replica": replica_id,
        "status": status,
        "lambda_single": None,
        "collective_coupling": None,
        "runtime_ps": None,
    }


def validate_replicas(
    output_dir: Path,
    verbose: bool = False,
) -> Tuple[int, int, int, List[Tuple[int, int]]]:
    """
    Validate that all expected replicas have completed.
    
    Returns:
        Tuple of (total_expected, total_success, total_failed, missing_list)
    """
    total_expected = 0
    total_success = 0
    total_failed = 0
    missing = []
    
    for n_molecules, expected_count in EXPECTED_REPLICAS.items():
        for replica_id in range(expected_count):
            total_expected += 1
            replica_dir = output_dir / f"N_{n_molecules}" / f"replica_{replica_id}"
            
            if not replica_dir.exists():
                missing.append((n_molecules, replica_id))
                continue
            
            status = get_replica_status(replica_dir)
            if status == "success":
                total_success += 1
            else:
                total_failed += 1
    
    return total_expected, total_success, total_failed, missing


def main() -> int:
    """Main entry point."""
    args = parse_args()
    
    output_dir = Path(args.output_dir)
    if not output_dir.exists():
        print(f"Error: Output directory does not exist: {output_dir}")
        return 1
    
    summary_file = output_dir / args.summary_file
    
    print("=" * 60)
    print("Finite-Size Scaling Study - Results Collection")
    print("=" * 60)
    print(f"Output directory: {output_dir}")
    print(f"Summary file: {summary_file}")
    print()
    
    # Validate if requested
    if args.validate:
        print("Validating replicas...")
        total_expected, total_success, total_failed, missing = validate_replicas(
            output_dir, args.verbose
        )
        
        print(f"  Expected: {total_expected}")
        print(f"  Success:  {total_success}")
        print(f"  Failed:   {total_failed}")
        print(f"  Missing:  {len(missing)}")
        
        if missing and args.verbose:
            print("\nMissing replicas:")
            for n, r in missing[:20]:  # Show first 20
                print(f"  N={n}, replica={r}")
            if len(missing) > 20:
                print(f"  ... and {len(missing) - 20} more")
        
        if total_success < total_expected:
            print(f"\nWARNING: Only {total_success}/{total_expected} replicas completed!")
            print("Consider re-running failed/missing replicas before collecting results.")
        print()
    
    # Collect all results
    print("Collecting results...")
    all_results = []
    
    for n_molecules in sorted(EXPECTED_REPLICAS.keys()):
        expected_count = EXPECTED_REPLICAS[n_molecules]
        success_count = 0
        
        if args.verbose:
            print(f"\nN={n_molecules} ({expected_count} replicas):")
        
        for replica_id in range(expected_count):
            data = collect_replica_data(
                output_dir, n_molecules, replica_id, args.verbose
            )
            if data:
                all_results.append(data)
                if data.get("status") == "success":
                    success_count += 1
        
        if not args.verbose:
            print(f"  N={n_molecules}: {success_count}/{expected_count} successful")
    
    # Write summary CSV
    print(f"\nWriting summary to: {summary_file}")
    
    fieldnames = [
        "N_molecules",
        "replica",
        "lambda_single",
        "collective_coupling",
        "runtime_ps",
        "device",
        "status",
    ]
    
    with summary_file.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in all_results:
            writer.writerow(row)
    
    # Summary statistics
    success_count = sum(1 for r in all_results if r.get("status") == "success")
    
    print()
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Total replicas collected: {len(all_results)}")
    print(f"Successful: {success_count}")
    print(f"Failed/Other: {len(all_results) - success_count}")
    print(f"Output file: {summary_file}")
    print("=" * 60)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
