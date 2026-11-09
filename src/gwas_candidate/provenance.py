"""Run provenance and audit-manifest helpers."""

import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy
import matplotlib
import yaml


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def executable_version(command):
    try:
        result = subprocess.run(
            [command, "-version"], text=True, capture_output=True, timeout=30
        )
    except (OSError, subprocess.TimeoutExpired):
        return "not_available"
    text = (result.stdout or result.stderr).strip().splitlines()
    return text[0] if text else "unknown"


def collect_versions():
    return {
        "python": platform.python_version(),
        "numpy": numpy.__version__,
        "matplotlib": matplotlib.__version__,
        "pyyaml": yaml.__version__,
        "blastp": executable_version("blastp"),
        "makeblastdb": executable_version("makeblastdb"),
    }


def _display_path(path, outdir):
    path = Path(path)
    try:
        return path.resolve().relative_to(Path(outdir).resolve()).as_posix()
    except ValueError:
        if path.is_absolute():
            return f"{path.parent.name}/{path.name}"
        return path.as_posix()


def _sanitize_config(value):
    if isinstance(value, dict):
        return {key: _sanitize_config(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_sanitize_config(item) for item in value]
    if isinstance(value, str):
        candidate = Path(value)
        if candidate.is_absolute():
            return f"{candidate.parent.name}/{candidate.name}"
    return value


def write_provenance(outdir, scenario, config, input_paths, commands, random_seed, reference_build):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    snapshot = outdir / "config_snapshot.yaml"
    snapshot.write_text(
        yaml.safe_dump(_sanitize_config(config), allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )
    manifest = {
        "scenario": scenario,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "random_seed": random_seed,
        "reference_genome_version": reference_build,
        "config_snapshot": snapshot.name,
        "input_sha256": {
            _display_path(path, outdir): sha256_file(path)
            for path in input_paths if Path(path).exists()
        },
        "software_versions": collect_versions(),
        "commands": commands,
    }
    path = outdir / "run_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return path
