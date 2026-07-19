#!/usr/bin/env python3
"""Prune aging_weak_lambda time-series data to a uniform sampling period."""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import h5py
import numpy as np

DEFAULT_TARGET_PERIOD_PS = 1.0
DEFAULT_ORIGIN_PS = 0.0
GRID_TOLERANCE_PS = 1.0e-3
VALUE_TOLERANCE = 1.0e-10

LAMBDA_DIRS = (
    "lambda0",
    "lambda0p01",
    "lambda0p016667",
    "lambda0p023333",
    "lambda0p03",
)


@dataclass(frozen=True)
class PruneStats:
    path: str
    layer: str
    n_in: int
    n_out: int
    size_in_bytes: int
    size_out_bytes: int
    target_period_ps: float
    elapsed_s: float
    status: str
    message: str = ""

    def to_dict(self) -> dict[str, object]:
        return {
            "path": self.path,
            "layer": self.layer,
            "n_in": self.n_in,
            "n_out": self.n_out,
            "size_in_bytes": self.size_in_bytes,
            "size_out_bytes": self.size_out_bytes,
            "target_period_ps": self.target_period_ps,
            "elapsed_s": self.elapsed_s,
            "status": self.status,
            "message": self.message,
        }


def trim_padded_trajectory(
    time_ps: np.ndarray,
    *series: np.ndarray,
    probe_index: int | None = None,
) -> tuple[np.ndarray, ...] | None:
    """Remove HDF5 trailing zero pads and leading all-zero samples."""
    time_arr = np.asarray(time_ps, dtype=float)
    arrays = [np.asarray(values, dtype=float) for values in series]
    if time_arr.size == 0 or not arrays:
        return None

    positive = time_arr > 0.0
    if not np.any(positive):
        return None

    end = int(np.where(positive)[0][-1]) + 1
    time_arr = time_arr[:end]
    arrays = [values[:end] for values in arrays]

    if np.any(np.diff(time_arr) < 0.0):
        order = np.argsort(time_arr)
        time_arr = time_arr[order]
        arrays = [values[order] for values in arrays]

    if probe_index is None:
        probe_index = 0
        for idx, values in enumerate(arrays):
            finite = values[np.isfinite(values)]
            if finite.size and np.any(finite != 0.0):
                probe_index = idx
                break

    leading = arrays[probe_index]
    valid = np.isfinite(leading) & (leading != 0.0)
    nonzero = np.where(valid)[0]
    if nonzero.size == 0:
        start = 0
    else:
        start = int(nonzero[0])
    trimmed = (time_arr[start:],) + tuple(values[start:] for values in arrays)
    return trimmed


def downsample_to_period(
    time_ps: np.ndarray,
    *series: np.ndarray,
    target_period_ps: float = DEFAULT_TARGET_PERIOD_PS,
    origin_ps: float = DEFAULT_ORIGIN_PS,
) -> tuple[np.ndarray, ...]:
    """Select samples nearest a uniform grid with spacing ``target_period_ps``."""
    time_arr = np.asarray(time_ps, dtype=float)
    arrays = [np.asarray(values, dtype=float) for values in series]
    if time_arr.size == 0:
        empty = (time_arr,) + tuple(values[:0] for values in arrays)
        return empty

    t_max = float(time_arr[-1])
    if t_max < origin_ps:
        empty = (time_arr[:0],) + tuple(values[:0] for values in arrays)
        return empty

    grid = np.arange(origin_ps, t_max + 0.5 * target_period_ps, target_period_ps)
    indices: list[int] = []
    grid_out: list[float] = []
    for target in grid:
        idx = int(np.argmin(np.abs(time_arr - target)))
        if indices and idx == indices[-1]:
            continue
        indices.append(idx)
        grid_out.append(float(target))

    idx_arr = np.asarray(indices, dtype=int)
    out_time = np.asarray(grid_out, dtype=float)
    downsampled = (out_time,) + tuple(values[idx_arr] for values in arrays)
    return downsampled


def is_uniform_grid(times: np.ndarray, period_ps: float, atol: float = GRID_TOLERANCE_PS) -> bool:
    """Return True when ``times`` is an evenly spaced 1-D grid."""
    times_arr = np.asarray(times, dtype=float)
    if times_arr.size < 2:
        return True
    diffs = np.diff(times_arr)
    return bool(np.allclose(diffs, period_ps, atol=atol))


def _iter_hdf5_datasets(group: h5py.Group, prefix: str = "") -> list[tuple[str, h5py.Dataset]]:
    datasets: list[tuple[str, h5py.Dataset]] = []
    for key, item in group.items():
        path = f"{prefix}/{key}" if prefix else key
        if isinstance(item, h5py.Dataset):
            datasets.append((path, item))
        elif isinstance(item, h5py.Group):
            datasets.extend(_iter_hdf5_datasets(item, path))
    return datasets


def _copy_dataset_attrs(src: h5py.Dataset, dst: h5py.Dataset) -> None:
    for attr_key, attr_val in src.attrs.items():
        dst.attrs[attr_key] = attr_val


def _write_group_tree_from_flat(
    src_group: h5py.Group,
    dst_group: h5py.Group,
    dataset_data: dict[str, np.ndarray],
) -> None:
    for key, item in src_group.items():
        if isinstance(item, h5py.Dataset):
            data = dataset_data[key]
            ds = dst_group.create_dataset(key, data=data, dtype=np.float64)
            _copy_dataset_attrs(item, ds)
        elif isinstance(item, h5py.Group):
            child = dst_group.create_group(key)
            nested = {
                path.split("/", 1)[1]: values
                for path, values in dataset_data.items()
                if path.startswith(f"{key}/")
            }
            _write_group_tree_from_flat(item, child, nested)


def prune_hdf5_observables(
    src: Path,
    dst: Path,
    target_period_ps: float = DEFAULT_TARGET_PERIOD_PS,
    origin_ps: float = DEFAULT_ORIGIN_PS,
) -> dict[str, object]:
    """Downsample all aligned observables in an HDF5 file."""
    with h5py.File(src, "r") as handle:
        if "time" not in handle:
            raise ValueError(f"Missing root time dataset: {src}")
        native_period = float(handle.attrs.get("output_period_ps", 0.01))
        time_raw = handle["time"][:]
        dataset_paths = _iter_hdf5_datasets(handle)
        path_order = [path for path, _ in dataset_paths if path != "time"]

        path_set = {path for path, _ in dataset_paths}
        probe_path = next(
            (
                candidate
                for candidate in (
                    "energies/system_total",
                    "energies/harmonic",
                    "energies/universe_total",
                )
                if candidate in path_set
            ),
            path_order[0],
        )
        trimmed = trim_padded_trajectory(time_raw, handle[probe_path][:])
        if trimmed is None:
            raise ValueError(f"No usable trajectory in {src}")
        trimmed_time = trimmed[0]
        start = int(np.where(time_raw == trimmed_time[0])[0][0])
        end = start + trimmed_time.size
        t_max = float(trimmed_time[-1])
        grid = np.arange(
            origin_ps,
            t_max + 0.5 * target_period_ps,
            target_period_ps,
        )
        index_list: list[int] = []
        for target in grid:
            rel = int(np.argmin(np.abs(trimmed_time - target)))
            idx = start + rel
            if index_list and idx == index_list[-1]:
                continue
            index_list.append(idx)
        index_arr = np.asarray(index_list, dtype=int)
        out_time = grid[: index_arr.size].astype(np.float64)
        out_map = {path: handle[path][index_arr] for path in path_order}
        out_map["time"] = out_time

    with h5py.File(src, "r") as src_handle, h5py.File(dst, "w", libver="latest") as dst_handle:
        for attr_key, attr_val in src_handle.attrs.items():
            if attr_key == "output_period_ps":
                dst_handle.attrs[attr_key] = target_period_ps
            else:
                dst_handle.attrs[attr_key] = attr_val
        if "output_period_ps" not in dst_handle.attrs:
            dst_handle.attrs["output_period_ps"] = target_period_ps

        for key, item in src_handle.items():
            if isinstance(item, h5py.Group):
                dst_group = dst_handle.create_group(key)
                nested = {
                    path.split("/", 1)[1]: out_map[path]
                    for path in out_map
                    if path.startswith(f"{key}/")
                }
                _write_group_tree_from_flat(item, dst_group, nested)
            elif isinstance(item, h5py.Dataset) and key == "time":
                ds = dst_handle.create_dataset("time", data=out_map["time"], dtype=np.float64)
                _copy_dataset_attrs(item, ds)

    return {
        "n_in": int(trimmed_time.size),
        "n_out": int(out_map["time"].shape[0]),
        "target_period_ps": target_period_ps,
        "t_min": float(out_map["time"][0]) if out_map["time"].size else None,
        "t_max": float(out_map["time"][-1]) if out_map["time"].size else None,
    }


def verify_pruned_hdf5(
    original: Path,
    pruned: Path,
    target_period_ps: float = DEFAULT_TARGET_PERIOD_PS,
    origin_ps: float = DEFAULT_ORIGIN_PS,
    n_spot_checks: int = 3,
) -> None:
    """Validate pruned HDF5 against stride-based expectations from the original."""
    with h5py.File(pruned, "r") as handle:
        if float(handle.attrs.get("output_period_ps", -1.0)) != target_period_ps:
            raise ValueError("output_period_ps mismatch")
        time_out = handle["time"][:]
        if time_out.size < 2:
            return
        if not is_uniform_grid(time_out, target_period_ps):
            raise ValueError("pruned time grid is not uniform at target period")

        for path, ds in _iter_hdf5_datasets(handle):
            if ds.shape[0] != time_out.shape[0]:
                raise ValueError(f"length mismatch for dataset {path}")

    with h5py.File(original, "r") as src_handle:
        native_period = float(src_handle.attrs.get("output_period_ps", 0.01))
        time_raw = src_handle["time"][:]
        paths = [path for path, _ in _iter_hdf5_datasets(src_handle) if path != "time"]
        path_set = set(paths)
        probe_path = next(
            (
                candidate
                for candidate in (
                    "energies/system_total",
                    "energies/harmonic",
                    "energies/universe_total",
                )
                if candidate in path_set
            ),
            paths[0],
        )
        trimmed = trim_padded_trajectory(time_raw, src_handle[probe_path][:])
        if trimmed is None:
            raise ValueError("original trajectory empty after trim")
        trimmed_time = trimmed[0]
        start = int(np.where(time_raw == trimmed_time[0])[0][0])
        t_max = float(trimmed_time[-1])
        grid = np.arange(
            origin_ps,
            t_max + 0.5 * target_period_ps,
            target_period_ps,
        )
        index_list: list[int] = []
        for target in grid:
            rel = int(np.argmin(np.abs(trimmed_time - target)))
            idx = start + rel
            if index_list and idx == index_list[-1]:
                continue
            index_list.append(idx)
        index_arr = np.asarray(index_list, dtype=int)
        check_paths = paths[: min(n_spot_checks, len(paths))]
        expected = [src_handle[path][index_arr] for path in check_paths]

    with h5py.File(pruned, "r") as pruned_handle:
        check_paths = paths[: min(n_spot_checks, len(paths))]
        for idx, path in enumerate(check_paths):
            expected_vals = expected[idx]
            actual_vals = pruned_handle[path][:]
            if not np.allclose(actual_vals, expected_vals, rtol=0.0, atol=VALUE_TOLERANCE):
                raise ValueError(f"spot-check failed for dataset {path}")


def read_fkt_series(path: Path) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Parse a per-replica F(k,t) text file into header, lag times, and values."""
    header: list[str] = []
    times: list[float] = []
    values: list[float] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            header.append(line.rstrip("\n"))
            continue
        parts = re.split(r"\s+", line.strip())
        if len(parts) < 2:
            continue
        try:
            lag = float(parts[0])
            value = float(parts[1])
        except ValueError:
            continue
        times.append(lag)
        values.append(value)
    return header, np.asarray(times, dtype=float), np.asarray(values, dtype=float)


def _dedupe_preserve_first(times: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if times.size == 0:
        return times, values
    order = np.argsort(times)
    times = times[order]
    values = values[order]
    keep: list[int] = []
    last: float | None = None
    for idx, lag in enumerate(times):
        if last is None or abs(lag - last) > GRID_TOLERANCE_PS:
            keep.append(idx)
            last = float(lag)
    idx_arr = np.asarray(keep, dtype=int)
    return times[idx_arr], values[idx_arr]


def prune_fkt_replica_file(
    src: Path,
    dst: Path,
    target_period_ps: float = DEFAULT_TARGET_PERIOD_PS,
    origin_ps: float = DEFAULT_ORIGIN_PS,
) -> dict[str, object]:
    """Resample a per-replica F(k,t) file onto a uniform lag grid."""
    header, times, values = read_fkt_series(src)
    if times.size == 0:
        raise ValueError(f"No F(k,t) samples in {src}")

    times, values = _dedupe_preserve_first(times, values)
    t_max = float(times[-1])
    grid = np.arange(origin_ps, t_max + 0.5 * target_period_ps, target_period_ps)
    out_values = np.interp(grid, times, values, left=np.nan, right=np.nan)
    valid = ~np.isnan(out_values)
    grid = grid[valid]
    out_values = out_values[valid]

    if not header:
        header = [
            "# F(k,t) correlation function",
            "# lag_time_ps\tF(k,t)",
        ]
    dst.parent.mkdir(parents=True, exist_ok=True)
    with dst.open("w", encoding="utf-8") as handle:
        for line in header:
            handle.write(f"{line}\n")
        for lag, value in zip(grid, out_values, strict=True):
            handle.write(f"{lag:.6f}\t{value:.8f}\n")

    return {
        "n_in": int(times.size),
        "n_out": int(grid.size),
        "target_period_ps": target_period_ps,
        "t_max": t_max,
    }


def _is_fkt_already_uniform(path: Path, target_period_ps: float) -> bool:
    try:
        _, times, _ = read_fkt_series(path)
    except OSError:
        return False
    if times.size < 2:
        return True
    return is_uniform_grid(times, target_period_ps)


def _replace_in_place(tmp_path: Path, final_path: Path) -> None:
    os.replace(tmp_path, final_path)


def prune_hdf5_file_in_place(
    path: Path,
    target_period_ps: float,
    dry_run: bool = False,
) -> PruneStats:
    start = time.perf_counter()
    size_in = path.stat().st_size
    tmp_path = path.with_suffix(path.suffix + ".prune_tmp")

    with h5py.File(path, "r") as handle:
        current_period = float(handle.attrs.get("output_period_ps", 0.01))
        n_in = int(np.count_nonzero(handle["time"][:] > 0.0))

    if abs(current_period - target_period_ps) < GRID_TOLERANCE_PS:
        return PruneStats(
            path=str(path),
            layer="hdf5",
            n_in=n_in,
            n_out=n_in,
            size_in_bytes=size_in,
            size_out_bytes=size_in,
            target_period_ps=target_period_ps,
            elapsed_s=time.perf_counter() - start,
            status="skipped",
            message="already at target period",
        )

    if dry_run:
        return PruneStats(
            path=str(path),
            layer="hdf5",
            n_in=n_in,
            n_out=max(1, int(n_in * current_period / target_period_ps)),
            size_in_bytes=size_in,
            size_out_bytes=int(size_in * current_period / target_period_ps),
            target_period_ps=target_period_ps,
            elapsed_s=time.perf_counter() - start,
            status="dry_run",
        )

    stats = prune_hdf5_observables(path, tmp_path, target_period_ps=target_period_ps)
    verify_pruned_hdf5(path, tmp_path, target_period_ps=target_period_ps)
    size_out = tmp_path.stat().st_size
    _replace_in_place(tmp_path, path)
    return PruneStats(
        path=str(path),
        layer="hdf5",
        n_in=int(stats["n_in"]),
        n_out=int(stats["n_out"]),
        size_in_bytes=size_in,
        size_out_bytes=size_out,
        target_period_ps=target_period_ps,
        elapsed_s=time.perf_counter() - start,
        status="ok",
    )


def prune_fkt_file_in_place(
    path: Path,
    target_period_ps: float,
    dry_run: bool = False,
) -> PruneStats:
    start = time.perf_counter()
    size_in = path.stat().st_size
    _, times, _ = read_fkt_series(path)
    n_in = int(times.size)

    if _is_fkt_already_uniform(path, target_period_ps):
        return PruneStats(
            path=str(path),
            layer="fkt",
            n_in=n_in,
            n_out=n_in,
            size_in_bytes=size_in,
            size_out_bytes=size_in,
            target_period_ps=target_period_ps,
            elapsed_s=time.perf_counter() - start,
            status="skipped",
            message="already uniform",
        )

    if dry_run:
        t_max = float(times.max()) if times.size else 0.0
        n_out = int(t_max / target_period_ps) + 1 if times.size else 0
        return PruneStats(
            path=str(path),
            layer="fkt",
            n_in=n_in,
            n_out=n_out,
            size_in_bytes=size_in,
            size_out_bytes=max(1, n_out * 32),
            target_period_ps=target_period_ps,
            elapsed_s=time.perf_counter() - start,
            status="dry_run",
        )

    tmp_path = path.with_suffix(path.suffix + ".prune_tmp")
    stats = prune_fkt_replica_file(path, tmp_path, target_period_ps=target_period_ps)
    _, out_times, _ = read_fkt_series(tmp_path)
    if out_times.size >= 2 and not is_uniform_grid(out_times, target_period_ps):
        tmp_path.unlink(missing_ok=True)
        raise ValueError(f"F(k,t) output grid not uniform: {path}")
    size_out = tmp_path.stat().st_size
    _replace_in_place(tmp_path, path)
    return PruneStats(
        path=str(path),
        layer="fkt",
        n_in=int(stats["n_in"]),
        n_out=int(stats["n_out"]),
        size_in_bytes=size_in,
        size_out_bytes=size_out,
        target_period_ps=target_period_ps,
        elapsed_s=time.perf_counter() - start,
        status="ok",
    )


def discover_hdf5_files(data_root: Path, lambda_id: str | None, replica: int | None) -> list[Path]:
    """Find observables HDF5 files under the campaign tree."""
    files: list[Path] = []
    search_roots: list[Path] = []
    if lambda_id:
        search_roots.append(data_root / lambda_id)
    else:
        search_roots.extend(data_root / name for name in LAMBDA_DIRS if (data_root / name).is_dir())
    search_roots.append(data_root / "analysis" / "pe_vs_T_calib")
    search_roots.append(data_root / "analysis" / "pe_vs_T_calib_rf")
    for lam_name in LAMBDA_DIRS:
        lam_root = data_root / lam_name
        if lam_root.is_dir():
            search_roots.append(lam_root)

    seen: set[Path] = set()
    for root in search_roots:
        if not root.is_dir():
            continue
        for path in sorted(root.rglob("observables_replica_*.h5")):
            if replica is not None:
                match = re.search(r"observables_replica_(\d+)\.h5$", path.name)
                if not match or int(match.group(1)) != replica:
                    continue
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            files.append(path)
    return sorted(files)


def discover_fkt_files(data_root: Path, lambda_id: str | None) -> list[Path]:
    files: list[Path] = []
    roots = [data_root / name for name in LAMBDA_DIRS if (data_root / name).is_dir()]
    if lambda_id:
        roots = [data_root / lambda_id]
    for root in roots:
        for run_dir in sorted(root.glob("cavity_coupling_*")):
            files.extend(sorted(run_dir.glob("prod-*_fkt_ref_*.txt")))
    return files


def append_manifest(manifest_path: Path, stats: PruneStats) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    payload = stats.to_dict()
    payload["timestamp"] = datetime.now(timezone.utc).isoformat()
    with manifest_path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload) + "\n")


def _prune_file_task(
    layer: str,
    path_str: str,
    target_period_ps: float,
    dry_run: bool,
) -> PruneStats:
    """Process one file in a worker process (replicas are independent)."""
    path = Path(path_str)
    if layer == "hdf5":
        return prune_hdf5_file_in_place(path, target_period_ps=target_period_ps, dry_run=dry_run)
    if layer == "fkt":
        return prune_fkt_file_in_place(path, target_period_ps=target_period_ps, dry_run=dry_run)
    raise ValueError(f"Unsupported layer: {layer}")


def run_layer(
    layer: str,
    data_root: Path,
    target_period_ps: float,
    dry_run: bool,
    lambda_id: str | None,
    replica: int | None,
    manifest_path: Path,
    max_files: int | None = None,
    max_workers: int = 1,
) -> list[PruneStats]:
    if layer == "hdf5":
        paths = discover_hdf5_files(data_root, lambda_id, replica)
        worker = prune_hdf5_file_in_place
    elif layer == "fkt":
        paths = discover_fkt_files(data_root, lambda_id)
        worker = prune_fkt_file_in_place
    else:
        raise ValueError(f"Unsupported layer: {layer}")

    if max_files is not None:
        paths = paths[:max_files]

    results: list[PruneStats] = []
    total = len(paths)
    if total == 0:
        return results

    workers = max(1, min(max_workers, total))
    if workers == 1:
        for idx, path in enumerate(paths, start=1):
            print(f"[{layer}] {idx}/{total} {path}", flush=True)
            try:
                stats = worker(path, target_period_ps=target_period_ps, dry_run=dry_run)
            except (FileNotFoundError, OSError) as exc:
                stats = PruneStats(
                    path=str(path),
                    layer=layer,
                    n_in=0,
                    n_out=0,
                    size_in_bytes=path.stat().st_size if path.exists() else 0,
                    size_out_bytes=0,
                    target_period_ps=target_period_ps,
                    elapsed_s=0.0,
                    status="error",
                    message=str(exc),
                )
            results.append(stats)
            append_manifest(manifest_path, stats)
        return results

    print(f"[{layer}] processing {total} files with {workers} workers", flush=True)
    completed = 0
    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_map = {
            executor.submit(
                _prune_file_task,
                layer,
                str(path),
                target_period_ps,
                dry_run,
            ): path
            for path in paths
        }
        for future in as_completed(future_map):
            path = future_map[future]
            completed += 1
            try:
                stats = future.result()
            except FileNotFoundError as exc:
                stats = PruneStats(
                    path=str(path),
                    layer=layer,
                    n_in=0,
                    n_out=0,
                    size_in_bytes=0,
                    size_out_bytes=0,
                    target_period_ps=target_period_ps,
                    elapsed_s=0.0,
                    status="error",
                    message=f"missing: {exc}",
                )
            except OSError as exc:
                stats = PruneStats(
                    path=str(path),
                    layer=layer,
                    n_in=0,
                    n_out=0,
                    size_in_bytes=path.stat().st_size if path.exists() else 0,
                    size_out_bytes=0,
                    target_period_ps=target_period_ps,
                    elapsed_s=0.0,
                    status="error",
                    message=str(exc),
                )
            except Exception as exc:
                stats = PruneStats(
                    path=str(path),
                    layer=layer,
                    n_in=0,
                    n_out=0,
                    size_in_bytes=path.stat().st_size if path.exists() else 0,
                    size_out_bytes=0,
                    target_period_ps=target_period_ps,
                    elapsed_s=0.0,
                    status="error",
                    message=f"{type(exc).__name__}: {exc}",
                )
            results.append(stats)
            append_manifest(manifest_path, stats)
            print(
                f"[{layer}] {completed}/{total} {stats.status} {path.name} "
                f"({stats.n_in}->{stats.n_out}) {stats.message}",
                flush=True,
            )
    return results


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Root of aging_weak_lambda data tree",
    )
    parser.add_argument(
        "--target-period-ps",
        type=float,
        default=DEFAULT_TARGET_PERIOD_PS,
        help="Target uniform sampling period in picoseconds",
    )
    parser.add_argument(
        "--layer",
        choices=("hdf5", "fkt", "all"),
        default="hdf5",
        help="Which data layer to prune",
    )
    parser.add_argument("--lambda-id", type=str, default=None, help="Restrict to one lambda directory")
    parser.add_argument("--replica", type=int, default=None, help="Restrict HDF5 to one replica index")
    parser.add_argument("--dry-run", action="store_true", help="Report actions without modifying files")
    parser.add_argument("--max-files", type=int, default=None, help="Process at most N files")
    parser.add_argument(
        "--max-workers",
        type=int,
        default=1,
        help="Parallel worker processes (replica files are independent)",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="JSONL manifest path (default: analysis/prune_1ps_manifest_<timestamp>.jsonl)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    data_root = args.data_root.resolve()
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest_path = args.manifest or (data_root / "analysis" / f"prune_1ps_manifest_{timestamp}.jsonl")

    layers: Iterable[str]
    if args.layer == "all":
        layers = ("hdf5", "fkt")
    else:
        layers = (args.layer,)

    for layer in layers:
        run_layer(
            layer=layer,
            data_root=data_root,
            target_period_ps=args.target_period_ps,
            dry_run=args.dry_run,
            lambda_id=args.lambda_id,
            replica=args.replica,
            manifest_path=manifest_path,
            max_files=args.max_files,
            max_workers=args.max_workers,
        )
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
