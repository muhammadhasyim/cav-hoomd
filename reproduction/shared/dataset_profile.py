"""Dataset profiles for reproduction scripts (paper archive vs in-repo aging data)."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Iterable, Optional

try:
    import yaml
except ImportError:  # pragma: no cover - optional dependency
    yaml = None  # type: ignore[assignment]

PROFILES_DIR = Path(__file__).resolve().parents[1] / "profiles"
REPO_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class CouplingEntry:
    """One coupling-strength run in a dataset profile."""

    id: str
    axis_value: float
    epsilon_tag: str
    run_dir: str
    analysis_tag: Optional[str] = None
    file_tag: Optional[str] = None

    @property
    def effective_tag(self) -> str:
        """Unique tag for staged dirs and time-series filenames."""
        return self.file_tag or self.epsilon_tag

    @property
    def staged_dir_name(self) -> str:
        return f"cavity_coupling_{self.effective_tag}_switch_200.0ps"


@dataclass(frozen=True)
class DatasetProfile:
    """Resolved paths and coupling catalogue for a reproduction dataset."""

    name: str
    data_root: Path
    staged_root: Path
    analysis_dir: Optional[Path]
    calibration: Optional[Path]
    relaxation_equilibrium: Optional[Path]
    figures_dir: Optional[Path]
    axis: str
    switch_time_ps: float
    couplings: tuple[CouplingEntry, ...]

    @property
    def time_series_dir(self) -> Path:
        return self.staged_root / "time_series_output"

    def coupling_by_id(self, coupling_id: str) -> CouplingEntry:
        for entry in self.couplings:
            if entry.id == coupling_id:
                return entry
        raise KeyError(f"Unknown coupling id {coupling_id!r} in profile {self.name!r}")

    def entry_for_axis_value(self, value: float, tol: float = 1e-9) -> CouplingEntry:
        for entry in self.couplings:
            if abs(entry.axis_value - value) <= tol:
                return entry
        raise KeyError(f"No coupling with axis value {value} in profile {self.name!r}")

    def tag_map(self) -> dict[float, str]:
        """Map axis value -> effective filename tag."""
        return {entry.axis_value: entry.effective_tag for entry in self.couplings}


def _resolve_path(value: Any, base: Path) -> Optional[Path]:
    if value is None or value == "null":
        return None
    path = Path(str(value))
    if not path.is_absolute():
        path = (base / path).resolve()
    return path


def _parse_coupling(raw: dict[str, Any], axis: str) -> CouplingEntry:
    if axis == "lambda":
        axis_value = float(raw["lambda"])
    elif axis == "epsilon":
        axis_value = float(raw["epsilon"])
    else:
        raise ValueError(f"Unsupported axis {axis!r}")

    return CouplingEntry(
        id=str(raw["id"]),
        axis_value=axis_value,
        epsilon_tag=str(raw["epsilon_tag"]),
        run_dir=str(raw["run_dir"]),
        analysis_tag=raw.get("analysis_tag"),
        file_tag=raw.get("file_tag"),
    )


def _parse_profile_dict(data: dict[str, Any]) -> DatasetProfile:
    axis = str(data["axis"])
    repo = REPO_ROOT
    data_root = _resolve_path(data["data_root"], repo)
    staged_root = _resolve_path(data.get("staged_root", data["data_root"]), repo)
    assert data_root is not None
    assert staged_root is not None

    couplings = tuple(_parse_coupling(item, axis) for item in data["couplings"])
    return DatasetProfile(
        name=str(data["name"]),
        data_root=data_root,
        staged_root=staged_root,
        analysis_dir=_resolve_path(data.get("analysis_dir"), repo),
        calibration=_resolve_path(data.get("calibration"), repo),
        relaxation_equilibrium=_resolve_path(data.get("relaxation_equilibrium"), repo),
        figures_dir=_resolve_path(data.get("figures_dir"), repo),
        axis=axis,
        switch_time_ps=float(data.get("switch_time_ps", 200.0)),
        couplings=couplings,
    )


def _parse_minimal_yaml(text: str) -> dict[str, Any]:
    """Parse the small YAML subset used by reproduction profile files."""
    root: dict[str, Any] = {}
    current_list: list[dict[str, Any]] | None = None
    current_item: dict[str, Any] | None = None
    list_key: str | None = None

    def _coerce(value: str) -> Any:
        value = value.strip()
        if value in {"null", "~", ""}:
            return None
        if value.lower() in {"true", "false"}:
            return value.lower() == "true"
        try:
            if any(ch in value for ch in ".eE"):
                return float(value)
            return int(value)
        except ValueError:
            return value

    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip(" "))
        stripped = line.strip()

        if stripped.startswith("- "):
            if current_list is None:
                raise ValueError(f"Unexpected list item: {stripped}")
            current_item = {}
            current_list.append(current_item)
            item_body = stripped[2:].strip()
            if item_body:
                key, value = item_body.split(":", 1)
                current_item[key.strip()] = _coerce(value)
            continue

        if ":" not in stripped:
            continue
        key, value = stripped.split(":", 1)
        key = key.strip()
        value = value.strip()

        if indent == 0:
            current_list = None
            current_item = None
            list_key = None
            if value == "":
                current_list = []
                root[key] = current_list
                list_key = key
            else:
                root[key] = _coerce(value)
            continue

        if indent >= 2 and current_item is not None and value != "":
            current_item[key] = _coerce(value)

    return root


def load_profile(name: str = "paper") -> DatasetProfile:
    """Load a dataset profile by name from ``reproduction/profiles/<name>.yaml``."""
    path = PROFILES_DIR / f"{name}.yaml"
    if not path.is_file():
        raise FileNotFoundError(f"Profile not found: {path}")

    if yaml is not None:
        with open(path, encoding="utf-8") as fh:
            data = yaml.safe_load(fh)
    else:
        data = _parse_minimal_yaml(path.read_text(encoding="utf-8"))
    profile = _parse_profile_dict(data)
    if profile.name != name:
        raise ValueError(f"Profile file {path} has name={profile.name!r}, expected {name!r}")
    return profile


def resolve_run_dir(profile: DatasetProfile, coupling_id: str) -> Path:
    """Absolute path to a coupling run directory under ``data_root``."""
    entry = profile.coupling_by_id(coupling_id)
    return profile.data_root / entry.run_dir


def staged_coupling_dir(profile: DatasetProfile, coupling_id: str) -> Path:
    """Absolute path to the flat staged coupling directory."""
    entry = profile.coupling_by_id(coupling_id)
    return profile.staged_root / entry.staged_dir_name


def coupling_values(profile: DatasetProfile) -> list[float]:
    """Return axis values (λ or ε) for all couplings in profile order."""
    return [entry.axis_value for entry in profile.couplings]


def coupling_tags(profile: DatasetProfile) -> dict[float, str]:
    """Return axis value -> effective filename tag mapping."""
    return profile.tag_map()


def epsilon_tag_for(profile: DatasetProfile, axis_value: float) -> str:
    """Filename tag for a given axis value."""
    return profile.entry_for_axis_value(axis_value).effective_tag


def staged_coupling_dirs(profile: DatasetProfile) -> list[Path]:
    """All flat staged coupling directory paths."""
    return [profile.staged_root / entry.staged_dir_name for entry in profile.couplings]


def add_profile_args(parser: argparse.ArgumentParser, default: str = "paper") -> None:
    """Register standard ``--profile`` and ``--staged-root`` CLI flags."""
    parser.add_argument(
        "--profile",
        default=default,
        help=f"Dataset profile name (default: {default})",
    )
    parser.add_argument(
        "--staged-root",
        type=Path,
        default=None,
        help="Override staged data root from the profile",
    )


def setup_profile(args: argparse.Namespace, default: str = "paper") -> DatasetProfile:
    """Load profile from parsed CLI args, optionally overriding staged_root."""
    profile_name = getattr(args, "profile", None) or default
    profile = load_profile(profile_name)
    staged_root = getattr(args, "staged_root", None)
    if staged_root is not None:
        profile = replace(profile, staged_root=Path(staged_root).resolve())
    return profile


def iter_coupling_entries(profile: DatasetProfile) -> Iterable[CouplingEntry]:
    """Iterate coupling entries in profile order."""
    yield from profile.couplings
