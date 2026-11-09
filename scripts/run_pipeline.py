#!/usr/bin/env python3
"""Run the V3 candidate-annotation pipeline from a YAML config."""

import argparse
import json
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gwas_candidate.pipeline import run_candidate_pipeline


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="YAML pipeline configuration")
    parser.add_argument("--workdir", required=True, help="Output directory")
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    base_dir = Path(config.get("base_dir", "."))
    if not base_dir.is_absolute():
        config["base_dir"] = str((config_path.parent / base_dir).resolve())
    result = run_candidate_pipeline(config, Path(args.workdir).resolve())
    print(json.dumps({
        "loci": len(result["loci"]),
        "candidate_evidence_rows": len(result["candidates"]),
        "report_html": str(result["report"]["html"]),
        "report_markdown": str(result["report"]["markdown"]),
        "manifest": str(result["manifest"]),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
