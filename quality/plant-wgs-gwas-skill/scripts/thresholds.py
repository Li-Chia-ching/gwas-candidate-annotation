#!/usr/bin/env python3
"""Single threshold resolver shared by GWAS plots and quality checks."""

from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass(frozen=True)
class ThresholdDecision:
    value: float
    method: str
    source: str
    meff: int | None = None


def _resolve_path(config_path, value):
    path = Path(value)
    if path.is_absolute():
        return path
    return Path(config_path).resolve().parent / path


def resolve_threshold(config_path=None, cli_threshold=None):
    """Resolve one threshold contract; an explicit CLI value wins and is recorded."""
    if cli_threshold is not None:
        if not (0 < cli_threshold <= 1):
            raise SystemExit("ERROR: threshold must be within (0, 1]")
        return ThresholdDecision(float(cli_threshold), "command_line", "--threshold")
    if not config_path:
        return ThresholdDecision(5e-8, "fixed_5e-8", "default")
    with open(config_path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    threshold_cfg = cfg.get("threshold", {})
    method = threshold_cfg.get("method", "fixed_5e-8")
    if method == "fixed_5e-8":
        value = float(threshold_cfg.get("value", 5e-8))
        return ThresholdDecision(value, method, "threshold.value")
    if method == "bonferroni_meff":
        source = threshold_cfg.get("meff_source")
        if source != "ld_pruned_marker_count":
            raise SystemExit(f"ERROR: unsupported meff_source: {source}")
        workdir = _resolve_path(config_path, cfg.get("paths", {}).get("workdir", ""))
        marker_file = workdir / "structure" / "pruned.prune.in"
        if not marker_file.exists():
            raise SystemExit(f"ERROR: M_eff marker file is missing: {marker_file}")
        markers = [line for line in marker_file.read_text(encoding="utf-8").splitlines() if line.strip()]
        if not markers:
            raise SystemExit(f"ERROR: M_eff marker file is empty: {marker_file}")
        meff = len(markers)
        return ThresholdDecision(0.05 / meff, method, str(marker_file), meff)
    raise SystemExit(f"ERROR: threshold method is not implemented in V3 P0: {method}")
