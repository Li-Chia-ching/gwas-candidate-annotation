"""CDS/protein extraction with strict structural QC."""

import gzip
from pathlib import Path

from .schemas import ProteinRecord


CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


class ProteinQCError(ValueError):
    pass


def open_any(path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else open(
        path, encoding="utf-8"
    )


def parse_gff_models(gff_path):
    genes = {}
    transcripts = {}
    cds = {}
    with open_any(gff_path) as handle:
        for line_no, line in enumerate(handle, 1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = dict(kv.split("=", 1) for kv in fields[8].split(";") if "=" in kv)
            feature = fields[2]
            if feature == "gene":
                gene_id = attrs.get("ID")
                if not gene_id:
                    raise ValueError(f"{gff_path}:{line_no}: gene lacks ID")
                genes[gene_id] = {
                    "chrom": fields[0], "start": int(fields[3]), "end": int(fields[4]),
                    "strand": fields[6], "transcripts": [],
                }
            elif feature in ("mRNA", "transcript"):
                transcript_id = attrs.get("ID")
                parent = attrs.get("Parent")
                if not transcript_id or not parent:
                    raise ValueError(f"{gff_path}:{line_no}: transcript requires ID and Parent")
                transcripts[transcript_id] = parent
                if parent in genes:
                    genes[parent]["transcripts"].append(transcript_id)
            elif feature == "CDS":
                parents = [x for x in attrs.get("Parent", "").split(",") if x]
                for parent in parents:
                    cds.setdefault(parent, []).append((int(fields[3]), int(fields[4])))
    return genes, transcripts, cds


def read_fasta(path):
    sequences = {}
    name = None
    chunks = []
    with open_any(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
        if name is not None:
            sequences[name] = "".join(chunks)
    return sequences


def revcomp(sequence):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def translate(sequence):
    sequence = sequence.upper()
    return "".join(CODON.get(sequence[i:i + 3], "X") for i in range(0, len(sequence) - 2, 3))


def qc_sequence(cds, protein):
    issues = []
    if len(cds) % 3 != 0:
        issues.append("cds_length_not_multiple_of_3")
    if cds[:3].upper() != "ATG":
        issues.append("missing_start_codon")
    if cds[-3:].upper() not in ("TAA", "TAG", "TGA"):
        issues.append("missing_stop_codon")
    if "*" in protein[:-1]:
        issues.append("internal_stop_codon")
    if "X" in protein:
        issues.append("unknown_amino_acid")
    return tuple(issues)


def extract_cds_protein(gff_path, genome_path, gene_ids, strict=True):
    genes, transcripts, cds = parse_gff_models(gff_path)
    genome = read_fasta(genome_path)
    records = []
    for gene_id in gene_ids:
        if gene_id not in genes:
            raise ValueError(f"gene {gene_id} is absent from GFF")
        gene = genes[gene_id]
        transcript_ids = [t for t in gene["transcripts"] if cds.get(t)]
        if not transcript_ids:
            raise ValueError(f"gene {gene_id} has no CDS-bearing transcript")
        transcript_id = max(
            transcript_ids,
            key=lambda t: sum(end - start + 1 for start, end in cds[t]),
        )
        chrom = gene["chrom"]
        if chrom not in genome:
            raise ValueError(f"chromosome {chrom} is absent from genome FASTA")
        sequence = genome[chrom]
        cds_seq = "".join(
            sequence[start - 1:end] for start, end in sorted(cds[transcript_id])
        )
        if gene["strand"] == "-":
            cds_seq = revcomp(cds_seq)
        protein = translate(cds_seq)
        issues = qc_sequence(cds_seq, protein)
        if strict and issues:
            raise ProteinQCError(f"{gene_id}: {', '.join(issues)}")
        records.append(ProteinRecord(
            gene_id=gene_id,
            transcript_id=transcript_id,
            chrom=chrom,
            start=gene["start"],
            end=gene["end"],
            strand=gene["strand"],
            cds_sequence=cds_seq,
            protein_sequence=protein,
            qc_status="OK" if not issues else "FAIL",
            qc_issues=issues,
        ))
    return records


def extract_proteome(gff_path, genome_path):
    genes, _, _ = parse_gff_models(gff_path)
    return extract_cds_protein(gff_path, genome_path, sorted(genes), strict=False)


def write_sequence_outputs(records, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cds_path = outdir / "candidate_cds.fasta"
    protein_path = outdir / "candidate_proteins.fasta"
    report_path = outdir / "extraction_report.tsv"
    with cds_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.gene_id}|{record.transcript_id}\n{record.cds_sequence}\n")
    with protein_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.gene_id}|{record.transcript_id}\n{record.protein_sequence}\n")
    with report_path.open("w", encoding="utf-8") as handle:
        handle.write("gene\ttranscript\tstrand\tcds_len\tprotein_len\tqc_status\tqc_issues\n")
        for record in records:
            handle.write(
                f"{record.gene_id}\t{record.transcript_id}\t{record.strand}\t"
                f"{len(record.cds_sequence)}\t{len(record.protein_sequence)}\t"
                f"{record.qc_status}\t{';'.join(record.qc_issues)}\n"
            )
    return cds_path, protein_path, report_path
