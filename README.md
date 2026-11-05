# gwas-candidate-annotation

An end-to-end Skill that takes GWAS significant signals to candidate-gene experimental prioritization. Designed for downstream analysis of crop resequencing GWAS; validated on an autotetraploid crop and equally applicable to diploids.

## What this Skill does

Input: GWAS results (rMVP/GAPIT-style signals.csv / summary statistics / ann.vcf). Output:

1. **Dual-model lead SNP list**: signals significant at the same position in two models (e.g., MLM + FarmCPU) → clustered into independent loci → one lead SNP per locus;
2. **Candidate gene table**: genes within a ±100 kb window of each lead SNP (from GFF);
3. **Candidate gene identity**: CDS extraction → protein translation → multi-database BLASTp (same species / close relatives / model plants / self-proteome) → domain and pathway annotation;
4. **Cloning feasibility assessment**: TE association detection (RVT/rve domains, fragmented hits), polyploid homolog copy inventory (prerequisite for copy-specific primers), and hormone-rescue gating (allowed only when the gene identity belongs to the corresponding hormone pathway).

## Quick start (10-minute walkthrough)

```bash
# 1. Dual-model lead screening (once per trait)
python scripts/screen_dual_model_leads.py \
    --signals-a trait.MLM_signals.csv --signals-b trait.FarmCPU_signals.csv \
    --model-a MLM --model-b FarmCPU \
    --gff reference.gff.gz --trait trait_name --outdir out/trait_name/

# 2. Extract candidate-gene CDS + proteins (--proteome builds a whole proteome for copy inventory)
python scripts/extract_cds_protein.py --gff reference.gff.gz \
    --genome genome.fasta.gz --genes gene1 gene2 --outdir out/anno/ --proteome

# 3. Multi-database BLASTp (requires local NCBI BLAST+; dbs.tsv format: database_name<TAB>protein_fasta_path)
python scripts/run_blastp.py --query out/anno/candidate_proteins.fasta \
    --dbs dbs.tsv --outdir out/anno/
python scripts/summarize_annotation.py --blastdir out/anno/blast_out \
    --outdir out/anno/ --self-db <self-proteome-db-name>

# 4. Write per-gene conclusions following the Stage 4 decision table in SKILL.md
```

## Directory structure

```
gwas-candidate-annotation/
├── SKILL.md                      Main document (YAML header + four-stage workflow + decision table)
├── README.md                     This file
├── scripts/
│   ├── screen_dual_model_leads.py   Dual-model lead SNP screening + window genes
│   ├── extract_cds_protein.py       CDS/protein extraction + whole-proteome construction
│   ├── run_blastp.py                Multi-database BLASTp (automatic database construction)
│   └── summarize_annotation.py      BLAST summary + homolog copy table
└── references/
    ├── workflow.md               Full workflow details + a validated plant-height trait case study (de-identified)
    ├── data_boundary.md          Interpretation boundaries and reporting language for truncated GWAS data
    └── pitfalls.md               Empirically verified pitfalls in polyploid cloning, TE detection, and tool environments
```

## Dependencies

- Python 3.9+ (standard library only)
- NCBI BLAST+ (blastp / makeblastdb; native Windows or Linux)
- Optional: HMMER (hmmscan) for Pfam domain scanning; when unavailable, inherit InterPro annotations from orthologs

## Three hard rules (see references/ for details)

1. **Coordinate system isolation**: use exactly one reference genome version per analysis; mixing versions is prohibited;
2. **Statistical significance ≠ clonability**: TE association (RVT/rve), fragmented hits, and zero-expression orthologs → all downgraded to KASP marker validation;
3. **Hormone-rescue gating**: a hormone-rescue experiment may be designed only when the annotated gene identity clearly belongs to that hormone pathway; defaulting by trait is prohibited ("plant height → GA" is a wrong default).

## Validation source

The workflow was distilled from a real downstream resequencing GWAS analysis in an autotetraploid crop (three quantitative traits × three association models, population size 200+, vendor-truncated delivery scenario). It includes complete annotation and tiered interpretation of two missense candidate genes (one highly conserved metabolic enzyme → clonable; one reverse transcriptase-related, likely TE-derived → downgraded). Specific gene identities and project data are not included in this Skill; see the internal records of the user's project for details.
