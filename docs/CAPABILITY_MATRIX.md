# Capability matrix

| Capability | Status | Evidence |
| --- | --- | --- |
| Dual-model summary filtering | Implemented | `test_p0_01_summary_filter.py` |
| Distance clusters without LD claims | Implemented | `test_p0_02_cluster_not_ld_locus.py` |
| Lead-lead LD merging of adjacent clusters | Implemented | `test_v31_ld_merge_rank.py`; opt-in via `ld_merge` + `vcf` |
| Locus ranking by nominal lead P | Implemented | `test_v31_ld_merge_rank.py`; `locus_rank`/`priority_peak` |
| Primary lead and cross-model support | Implemented | `test_p0_03_primary_lead.py` |
| REF/ALT, effect, reference-build gates | Implemented | `test_p0_04_allele_effect_reference.py` |
| PLINK2/GEMMA field adapters | Implemented | `test_p0_05_native_adapters.py` |
| GEMMA input preparation | Constructed-input tested | `test_p0_06_gemma_inputs.py`; real-tool validation pending |
| Shared M_eff threshold | Implemented | `test_p0_07_meff_threshold.py` |
| Full/truncated mode | Implemented | `test_p0_08_data_mode.py` |
| GFF window/direct/gene-body candidates | Implemented | `test_migration_candidates.py` |
| Strict CDS/protein QC | Implemented | `test_migration_sequence.py` |
| Multi-database BLASTp | Implemented | `test_migration_homology.py`; E2E uses BLAST+ 2.17.0 |
| Self-BLAST homoeolog candidates | Implemented | `test_06_polyploid_self_blast_end_to_end` |
| TE and hormone gates | Implemented | `test_migration_gates.py` |
| HTML/Markdown/TSV/PDF report | Implemented | `test_report_generator.py` |
| Statistical fine-mapping | Planned | Not implemented |
| Expression/eQTL/haplotype adapters | Planned | Not implemented |
| Pfam/InterPro scanning | Planned | User-supplied `domains.tsv` required |
| True ortholog inference | Planned | BLAST hits remain homology candidates |
| CNV/true copy number | Planned | Self-BLAST hits remain homoeolog candidates |
| MLMM/FarmCPU/BLINK execution | External dependency | Not implemented in repository |
| `ann.vcf` ANN/HGVS parsing | Planned | Not implemented |
| PLINK2/GEMMA real-data execution | External dependency, not CI-covered | Constructed-input test only |
