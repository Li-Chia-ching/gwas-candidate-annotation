# Changelog

## v3.1.0 - 2026-11-09

### Added

- `merge_clusters_by_ld`: single-linkage merging of adjacent distance clusters whose lead-SNP dosage r2 exceeds a configurable threshold (default 0.5) from a biallelic VCF. The merged locus takes the minimum-P variant as nominal lead and is marked `cluster_type=ld_merged_cluster`, `ld_locus_status=merged_by_ld`, `merged_by_ld`, `merged_cluster_count`.
- `rank_loci`: loci ranked by nominal lead P (ascending) with `locus_rank` and `priority_peak` (rank 1); ties break deterministically on chromosome/position.
- Pipeline config keys `ld_merge` (opt-in; hard-fails without `vcf`) and `ld_merge_r2` (default 0.5). Ranking is applied in every run.
- `tables/loci.tsv` gains `locus_rank`, `priority_peak`, `lead_p`, `merged_by_ld`, `merged_cluster_count`; the report gains a ranked-loci table and an explicit "nominal lead; conditional analysis not implemented" note in Markdown and HTML.
- 12 new regression tests: LD merge, low-LD no-merge, no-VCF degradation, missing lead genotypes, transitive merge with lead recomputation, chromosome naming fallback, r2 threshold validation, ranking order and tie-breaks, pipeline integration, and default no-merge behaviour.

### Notes

- Lead SNPs remain nominal leads (minimum P within a locus); conditional analysis (COJO/stepwise) for multiple independent signals within one locus is not implemented and is labelled as such in every report.
- Loci whose lead genotypes are absent from the VCF keep `ld_locus_status=not_evaluated` and are never merged silently.
- Distance clustering without a VCF is unchanged from v3.0.0; P0 behaviour and the legacy screening script are untouched.

### Validation

- 44/44 tests passed: 14 P0 regression, 10 migration regression, 2 report, 6 constructed-data end-to-end, 12 v3.1 LD-merge/rank regression.

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
