# Release manifest

Version: 3.1.0

## Included

- Root metadata: `LICENSE`, `VERSION`, `CHANGELOG.md`, `THIRD_PARTY.md`, `README.md`, `SKILL.md`, `pyproject.toml`, `requirements.txt`, `Makefile`, `.gitignore`
- CLI: `scripts/run_pipeline.py`
- V3 package: `src/gwas_candidate/`
- Tests: `tests/`
- Constructed quickstart data: `examples/minimal/`
- Documentation: `docs/`
- Audited historical source: `legacy/gwas-candidate-annotation/`
- Audited quality framework: `quality/plant-wgs-gwas-skill/`

## Excluded from the ZIP

- `.git/`
- `artifacts/`
- `out/`
- `__pycache__/`, `.pytest_cache/`, `.venv/`
- `.idea/`, `.vscode/`, `.DS_Store`
- Local credentials, keys, environment files, and personal paths

The exact ZIP entry count and SHA-256 are recorded in the release notes produced after packaging.
