# Changelog

## v3.0.0 - 2026-11-09

### Added

- V3 adapter package for dual-model association inputs, GFF candidates, CDS/protein QC, BLAST homology candidates, experimental gates, LD heatmaps, and evidence reports.
- `scripts/run_pipeline.py` CLI and `examples/minimal/` constructed-data quickstart.
- `report.html`, `report.md`, standardized TSV tables, and independent 600-dpi PDF figures.
- Run manifests with relative input paths, SHA-256 checksums, config snapshots, software versions, commands, random seed, and reference genome version.

### Fixed

- Full summary statistics are filtered before dual-model intersection.
- Distance clusters are no longer labeled independent LD loci.
- Primary-model lead SNP is separated from cross-model support.
- REF/ALT, multiallelic conflicts, effect direction, and reference-build mismatches block execution.
- Native PLINK2 and GEMMA result columns are adapted for diagnostics.
- GEMMA phenotype, numeric covariate, and kinship inputs are sample-aligned.
- Plots and checks share the same Bonferroni M_eff threshold resolver.
- Full and truncated data modes are explicit; truncated inputs disable Q-Q and lambda.
- Report and manifest paths are relative to avoid local-path leakage.

### Validation

- 14 P0 regression tests.
- 10 migration regression tests.
- 2 report tests.
- 6 constructed-data end-to-end scenarios.
- 32/32 total tests passed with Python 3.13.14, NumPy 2.5.4, Matplotlib 3.11.2, PyYAML 6.0.3, and NCBI BLAST+ 2.17.0.
