# Minimal constructed-data example

All records are synthetic and safe to redistribute.

Run from the repository root:

```bash
python scripts/run_pipeline.py --config examples/minimal/config.yaml --workdir out/minimal
```

Expected outputs:

- `out/minimal/report/report.html`
- `out/minimal/report/report.md`
- `out/minimal/report/tables/*.tsv`
- `out/minimal/report/figures/*.pdf`
- `out/minimal/run_manifest.json`
- `out/minimal/config_snapshot.yaml`

Expected result summary: one distance cluster, one candidate gene, one BLAST homology candidate, Priority B, Evidence L2.
