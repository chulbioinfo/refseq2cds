#!/usr/bin/env python3
"""
Build per-gene human CDS-position to genomic-position matrices.

Input assumptions:
- selection/representative_cds.normalized.parquet exists
- fastas/manifest.tsv exists
- raw/assembly_packages/human contains the assembly-exact GFF3 and sequence report

Output:
- human_cds_matrices/{SYMBOL}.human_cds_genomic_matrix.tsv.gz
- human_cds_matrices/manifest.tsv
- reports/human_cds_position_matrices.summary.json

Each matrix is nucleotide-level. Rows are ordered in the human CDS orientation
used in fastas/{SYMBOL}.fasta after terminal stop normalization.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.parse import unquote

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
ROOT = PROJECT_ROOT
OUTPUT_ROOT = PROJECT_ROOT
MANIFEST = PROJECT_ROOT / "config" / "species_manifest.tsv"
RAW_ASSEMBLIES = PROJECT_ROOT / "raw" / "assembly_packages"
SELECTION = PROJECT_ROOT / "selection" / "representative_cds.normalized.parquet"
FASTA_MANIFEST = PROJECT_ROOT / "fastas" / "manifest.tsv"
OUTDIR = PROJECT_ROOT / "human_cds_matrices"
REPORTS = PROJECT_ROOT / "reports"

GENEID_RE = re.compile(r"GeneID:(\d+)")
TRANSCRIPT_RE = re.compile(r"\b[NUX][MR]_\d+(?:\.\d+)?\b")
PROTEIN_RE = re.compile(r"\b[NX]P_\d+(?:\.\d+)?\b")


@dataclass
class CdsFeature:
    seqid: str
    start: int
    end: int
    strand: str
    phase: str
    geneid: int
    transcript_accession: str
    protein_accession: str


def log(message: str) -> None:
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def configure_paths(
    *,
    manifest: Optional[str] = None,
    input_root: Optional[str] = None,
    output_root: Optional[str] = None,
) -> None:
    global ROOT, OUTPUT_ROOT, MANIFEST, RAW_ASSEMBLIES, SELECTION, FASTA_MANIFEST, OUTDIR, REPORTS
    OUTPUT_ROOT = Path(output_root).resolve() if output_root else PROJECT_ROOT
    ROOT = OUTPUT_ROOT
    raw_root = Path(input_root).resolve() if input_root else OUTPUT_ROOT / "raw"
    MANIFEST = Path(manifest).resolve() if manifest else PROJECT_ROOT / "config" / "species_manifest.tsv"
    RAW_ASSEMBLIES = raw_root / "assembly_packages"
    SELECTION = OUTPUT_ROOT / "selection" / "representative_cds.normalized.parquet"
    FASTA_MANIFEST = OUTPUT_ROOT / "fastas" / "manifest.tsv"
    OUTDIR = OUTPUT_ROOT / "human_cds_matrices"
    REPORTS = OUTPUT_ROOT / "reports"


def rel_path(path: Path) -> str:
    for base in [OUTPUT_ROOT, PROJECT_ROOT]:
        try:
            return str(path.relative_to(base))
        except ValueError:
            continue
    return str(path)


def locate_file(package_dir: Path, predicate) -> Path:
    candidates = sorted([p for p in package_dir.rglob("*") if p.is_file() and predicate(p)])
    if not candidates:
        raise FileNotFoundError(f"No matching file under {package_dir}")
    return candidates[0]


def parse_attrs(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_text.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        else:
            key, value = part, ""
        attrs[unquote(key)] = unquote(value)
    return attrs


def accession_from_text(text: str, regex: re.Pattern[str]) -> str:
    match = regex.search(text or "")
    return match.group(0) if match else ""


def extract_geneid(attrs: Dict[str, str]) -> Optional[int]:
    for key in ("Dbxref", "db_xref"):
        match = GENEID_RE.search(attrs.get(key, ""))
        if match:
            return int(match.group(1))
    return None


def parse_sequence_report(path: Path) -> Dict[str, Dict[str, str]]:
    seqid_to_info: Dict[str, Dict[str, str]] = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            rec = json.loads(line)
            refseq = rec.get("refseqAccession") or ""
            if not refseq:
                continue
            seqid_to_info[refseq] = {
                "ucsc_chrom": rec.get("ucscStyleName") or refseq,
                "chr_name": rec.get("chrName") or rec.get("sequenceName") or "",
                "assembly_unit": rec.get("assemblyUnit") or "",
                "role": rec.get("role") or "",
            }
    return seqid_to_info


def parse_assembly_report(path: Path) -> Dict[str, str]:
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            rec = json.loads(line)
            assembly_info = rec.get("assemblyInfo") or rec.get("assembly_info") or {}
            annotation_info = rec.get("annotationInfo") or rec.get("annotation_info") or {}
            return {
                "gcf_accession": rec.get("accession") or "GCF_009914755.1",
                "assembly_name": assembly_info.get("assemblyName") or assembly_info.get("assembly_name") or "",
                "annotation_release_id": (
                    annotation_info.get("name")
                    or annotation_info.get("release_id")
                    or annotation_info.get("releaseId")
                    or ""
                ),
                "annotation_release_date": (
                    annotation_info.get("releaseDate")
                    or annotation_info.get("release_date")
                    or ""
                ),
            }
    return {
        "gcf_accession": "GCF_009914755.1",
        "assembly_name": "",
        "annotation_release_id": "",
        "annotation_release_date": "",
    }


def parse_human_cds_features(gff_path: Path) -> Dict[Tuple[int, str, str], List[CdsFeature]]:
    by_key: Dict[Tuple[int, str, str], List[CdsFeature]] = defaultdict(list)
    with gff_path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            seqid, _source, _ftype, start, end, _score, strand, phase, attr_text = parts
            attrs = parse_attrs(attr_text)
            geneid = extract_geneid(attrs)
            if geneid is None:
                continue
            parent = attrs.get("Parent", "")
            transcript = (
                accession_from_text(parent, TRANSCRIPT_RE)
                or accession_from_text(attrs.get("transcript_id", ""), TRANSCRIPT_RE)
                or accession_from_text(attrs.get("Dbxref", ""), TRANSCRIPT_RE)
            )
            protein = (
                attrs.get("protein_id", "")
                or accession_from_text(attrs.get("Name", ""), PROTEIN_RE)
                or accession_from_text(attrs.get("Dbxref", ""), PROTEIN_RE)
            )
            if not protein:
                continue
            feature = CdsFeature(
                seqid=seqid,
                start=int(start),
                end=int(end),
                strand=strand,
                phase=phase,
                geneid=geneid,
                transcript_accession=transcript,
                protein_accession=protein,
            )
            by_key[(geneid, transcript, protein)].append(feature)
    return by_key


def ordered_features(features: List[CdsFeature]) -> List[CdsFeature]:
    if not features:
        return []
    strand = features[0].strand
    if strand == "-":
        return sorted(features, key=lambda f: (f.seqid, f.start, f.end), reverse=True)
    return sorted(features, key=lambda f: (f.seqid, f.start, f.end))


def expand_genomic_positions(features: List[CdsFeature]) -> List[Tuple[int, CdsFeature, int]]:
    """Return transcript-order (genomic_position, feature, exon_rank) tuples."""
    ordered = ordered_features(features)
    positions: List[Tuple[int, CdsFeature, int]] = []
    for exon_rank, feat in enumerate(ordered, start=1):
        if feat.strand == "-":
            iterator = range(feat.end, feat.start - 1, -1)
        else:
            iterator = range(feat.start, feat.end + 1)
        for pos in iterator:
            positions.append((pos, feat, exon_rank))
    return positions


def safe_output_stem(fasta_path: str, human_symbol: str, human_geneid: int) -> str:
    if fasta_path:
        stem = Path(fasta_path).stem
        if stem:
            return stem
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(human_symbol)).strip("._")
    return cleaned or f"GeneID_{human_geneid}"


def selected_human_symbol(selected: Dict[str, Any]) -> str:
    return str(
        selected.get("species_symbol")
        or selected.get("actual_human_symbol")
        or selected.get("human_symbol")
        or f"GeneID_{selected.get('GeneID')}"
    )


def write_matrix(
    out_path: Path,
    selected: Dict[str, Any],
    positions: List[Tuple[int, CdsFeature, int]],
    seqid_info: Dict[str, Dict[str, str]],
    assembly_meta: Dict[str, str],
) -> Dict[str, Any]:
    normalized_seq = str(selected["normalized_sequence"]).upper()
    human_symbol = selected_human_symbol(selected)
    human_geneid = int(selected["GeneID"])
    terminal_stop_removed = bool(selected.get("terminal_stop_removed", False))
    original_position_count = len(positions)
    expected_length = int(selected["cds_length_normalized"])
    if len(normalized_seq) != expected_length:
        raise ValueError(
            f"{selected['family_id']} normalized sequence length mismatch: "
            f"sequence={len(normalized_seq)} metadata={expected_length}"
        )
    if original_position_count < expected_length:
        raise ValueError(
            f"{selected['family_id']} has fewer genomic CDS positions ({original_position_count}) "
            f"than normalized CDS length ({expected_length})"
        )
    positions = positions[:expected_length]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "family_id",
        "human_symbol",
        "human_GeneID",
        "transcript_accession",
        "protein_accession",
        "cds_nt_pos_1based",
        "codon_index_1based",
        "codon_offset_1based",
        "cds_base_transcript_orientation",
        "refseq_accession",
        "ucsc_chrom",
        "genomic_pos_1based",
        "bed_start_0based",
        "bed_end_0based",
        "strand",
        "cds_exon_rank",
        "cds_exon_start_1based",
        "cds_exon_end_1based",
        "gff_phase",
        "terminal_stop_removed",
        "gcf_accession",
        "assembly_name",
        "annotation_release_id",
        "annotation_release_date",
    ]
    with gzip.open(out_path, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for i, ((genomic_pos, feat, exon_rank), base) in enumerate(
            zip(positions, normalized_seq), start=1
        ):
            info = seqid_info.get(feat.seqid, {})
            writer.writerow(
                {
                    "family_id": selected["family_id"],
                    "human_symbol": human_symbol,
                    "human_GeneID": human_geneid,
                    "transcript_accession": selected["transcript_accession"],
                    "protein_accession": selected["protein_accession"],
                    "cds_nt_pos_1based": i,
                    "codon_index_1based": ((i - 1) // 3) + 1,
                    "codon_offset_1based": ((i - 1) % 3) + 1,
                    "cds_base_transcript_orientation": base,
                    "refseq_accession": feat.seqid,
                    "ucsc_chrom": info.get("ucsc_chrom", feat.seqid),
                    "genomic_pos_1based": genomic_pos,
                    "bed_start_0based": genomic_pos - 1,
                    "bed_end_0based": genomic_pos,
                    "strand": feat.strand,
                    "cds_exon_rank": exon_rank,
                    "cds_exon_start_1based": feat.start,
                    "cds_exon_end_1based": feat.end,
                    "gff_phase": feat.phase,
                    "terminal_stop_removed": terminal_stop_removed,
                    "gcf_accession": assembly_meta["gcf_accession"],
                    "assembly_name": assembly_meta["assembly_name"],
                    "annotation_release_id": assembly_meta["annotation_release_id"],
                    "annotation_release_date": assembly_meta["annotation_release_date"],
                }
            )
    return {
        "matrix_rows": expected_length,
        "original_cds_genomic_positions": original_position_count,
        "terminal_stop_removed": terminal_stop_removed,
        "strand": features_strand(positions),
    }


def features_strand(positions: List[Tuple[int, CdsFeature, int]]) -> str:
    strands = {feat.strand for _pos, feat, _rank in positions}
    return ",".join(sorted(strands))


def read_human_manifest_row() -> Dict[str, str]:
    with MANIFEST.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if int(row["taxid"]) == 9606:
                return row
    raise RuntimeError("config/species_manifest.tsv must include human taxid 9606 to build human CDS matrices")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default=str(MANIFEST), help="Species manifest TSV")
    parser.add_argument("--input-root", help="Input root containing assembly_packages/")
    parser.add_argument("--output-root", help="Output root containing selection/ and fastas/")
    parser.add_argument("--force", action="store_true", help="Rebuild existing matrix files")
    args = parser.parse_args()

    configure_paths(manifest=args.manifest, input_root=args.input_root, output_root=args.output_root)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)

    human_manifest = read_human_manifest_row()
    human_token = human_manifest["token"]
    human_package = RAW_ASSEMBLIES / human_token

    gff_path = locate_file(human_package, lambda p: p.name == "genomic.gff" or p.suffix == ".gff3")
    seq_report_path = locate_file(human_package, lambda p: p.name == "sequence_report.jsonl")
    assembly_report_path = locate_file(human_package, lambda p: p.name == "assembly_data_report.jsonl")

    log("Parsing human sequence report")
    seqid_info = parse_sequence_report(seq_report_path)
    assembly_meta = parse_assembly_report(assembly_report_path)

    log("Parsing human CDS features from GFF3")
    cds_features = parse_human_cds_features(gff_path)

    selected = pd.read_parquet(SELECTION)
    selected = selected[
        (selected["token"] == human_token)
        & (selected["family_cds_qc_passed"] == True)
    ].copy()
    fasta_manifest = pd.read_csv(FASTA_MANIFEST, sep="\t")
    fasta_by_family = fasta_manifest.set_index("family_id")["fasta_path"].to_dict()
    selected = selected[selected["family_id"].isin(fasta_by_family)].copy()

    manifest_rows: List[Dict[str, Any]] = []
    failures: List[Dict[str, Any]] = []
    total_rows = 0

    for idx, row in enumerate(selected.to_dict("records"), start=1):
        if idx == 1 or idx % 1000 == 0:
            log(f"Writing matrices: {idx:,}/{len(selected):,}")
        family_id = row["family_id"]
        key = (
            int(row["GeneID"]),
            str(row.get("transcript_accession", "")),
            str(row.get("protein_accession", "")),
        )
        features = cds_features.get(key)
        if not features:
            failures.append(
                {
                    "family_id": family_id,
                    "human_symbol": selected_human_symbol(row),
                    "human_GeneID": row.get("GeneID", ""),
                    "transcript_accession": row.get("transcript_accession", ""),
                    "protein_accession": row.get("protein_accession", ""),
                    "reason": "no_matching_gff_cds_features",
                }
            )
            continue
        positions = expand_genomic_positions(features)
        stem = safe_output_stem(fasta_by_family.get(family_id, ""), selected_human_symbol(row), int(row["GeneID"]))
        out_path = OUTDIR / f"{stem}.human_cds_genomic_matrix.tsv.gz"
        if out_path.exists() and not args.force:
            matrix_rows = int(row["cds_length_normalized"])
            stats = {
                "matrix_rows": matrix_rows,
                "original_cds_genomic_positions": len(positions),
                "terminal_stop_removed": bool(row.get("terminal_stop_removed", False)),
                "strand": features[0].strand,
            }
        else:
            try:
                stats = write_matrix(out_path, row, positions, seqid_info, assembly_meta)
            except Exception as exc:
                failures.append(
                    {
                        "family_id": family_id,
                        "human_symbol": selected_human_symbol(row),
                        "human_GeneID": row.get("GeneID", ""),
                        "transcript_accession": row.get("transcript_accession", ""),
                        "protein_accession": row.get("protein_accession", ""),
                        "reason": f"write_failed:{exc}",
                    }
                )
                continue
        total_rows += int(stats["matrix_rows"])
        manifest_rows.append(
            {
                "family_id": family_id,
                "human_symbol": selected_human_symbol(row),
                "human_GeneID": int(row["GeneID"]),
                "reference_symbol": row.get("reference_symbol", row.get("human_symbol", "")),
                "reference_GeneID": row.get("reference_GeneID", row.get("human_GeneID", "")),
                "reference_taxid": row.get("reference_taxid", ""),
                "reference_token": row.get("reference_token", ""),
                "transcript_accession": row.get("transcript_accession", ""),
                "protein_accession": row.get("protein_accession", ""),
                "matrix_path": rel_path(out_path),
                "matrix_rows": int(stats["matrix_rows"]),
                "original_cds_genomic_positions": int(stats["original_cds_genomic_positions"]),
                "terminal_stop_removed": bool(stats["terminal_stop_removed"]),
                "strand": stats["strand"],
                "gcf_accession": assembly_meta["gcf_accession"],
                "annotation_release_id": assembly_meta["annotation_release_id"],
            }
        )

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_df.to_csv(OUTDIR / "manifest.tsv", sep="\t", index=False)
    failures_df = pd.DataFrame(failures)
    failures_df.to_csv(OUTDIR / "failed.tsv", sep="\t", index=False)

    summary = {
        "generated_at": datetime.now().isoformat(),
        "input_human_families": int(len(selected)),
        "matrix_files": int(len(manifest_df)),
        "failed": int(len(failures_df)),
        "total_matrix_rows": int(total_rows),
        "output_dir": rel_path(OUTDIR),
        "gff_path": rel_path(gff_path),
        "sequence_report_path": rel_path(seq_report_path),
        **assembly_meta,
    }
    with (REPORTS / "human_cds_position_matrices.summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
        fh.write("\n")
    log(
        f"Done: matrix_files={summary['matrix_files']:,}; failed={summary['failed']:,}; "
        f"rows={summary['total_matrix_rows']:,}"
    )


if __name__ == "__main__":
    main()
