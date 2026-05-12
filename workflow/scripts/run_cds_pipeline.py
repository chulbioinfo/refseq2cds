#!/usr/bin/env python3
"""
Assembly-exact NCBI singleton CDS FASTA builder.

This driver implements the CDS-generation half of
NCBI_ortholog.plan_20260511.md:

1. freeze NCBI Gene FTP bulk files
2. download frozen GCF genome annotation packages
3. index assembly-exact GFF3/CDS/protein/RNA data
4. build ortholog components from gene_orthologs.gz
5. select one representative CDS per family/species
6. QC and write fastas/{HUMAN_SYMBOL}.fasta

The script is intentionally restartable: existing completed stage outputs are
reused unless --force is supplied.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple
from urllib.parse import unquote

import networkx as nx
import pandas as pd
from Bio import SeqIO
from Bio.Data.CodonTable import unambiguous_dna_by_id


PROJECT_ROOT = Path(__file__).resolve().parents[2]
ROOT = PROJECT_ROOT
OUTPUT_ROOT = PROJECT_ROOT
MANIFEST = PROJECT_ROOT / "config" / "species_manifest.tsv"
DATASETS = PROJECT_ROOT / "bin" / "datasets"
RAW_BULK = PROJECT_ROOT / "raw" / "ncbi_bulk"
RAW_ASSEMBLIES = PROJECT_ROOT / "raw" / "assembly_packages"
INDEXES = PROJECT_ROOT / "indexes"
ORTHOLOGY = PROJECT_ROOT / "orthology"
SELECTION = PROJECT_ROOT / "selection"
QC = PROJECT_ROOT / "qc"
FASTAS = PROJECT_ROOT / "fastas"
REPORTS = PROJECT_ROOT / "reports"
LOGS = PROJECT_ROOT / "logs"
RESULTS = PROJECT_ROOT / "results"
GENEID_FILTER_BY_TAXID: Optional[Dict[int, set[int]]] = None

STANDARD_TABLE = unambiguous_dna_by_id[1]
STOP_CODONS = set(STANDARD_TABLE.stop_codons)
DNA_RE = re.compile(r"^[ACGTNacgtn]+$")
TRANSCRIPT_RE = re.compile(r"\b[NUX][MR]_\d+(?:\.\d+)?\b")
PROTEIN_RE = re.compile(r"\b[NX]P_\d+(?:\.\d+)?\b")
GENEID_RE = re.compile(r"GeneID:(\d+)")


@dataclass(frozen=True)
class Species:
    token: str
    taxid: int
    scientific_name: str
    common_name: str
    gcf_accession: str
    outgroup: bool


def log(message: str) -> None:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {message}", flush=True)


def configure_paths(
    *,
    manifest: Optional[str] = None,
    input_root: Optional[str] = None,
    output_root: Optional[str] = None,
    datasets: Optional[str] = None,
) -> None:
    """Configure input/output roots without breaking the repository defaults."""
    global ROOT, OUTPUT_ROOT, MANIFEST, DATASETS
    global RAW_BULK, RAW_ASSEMBLIES, INDEXES, ORTHOLOGY, SELECTION, QC, FASTAS, REPORTS, LOGS, RESULTS

    OUTPUT_ROOT = Path(output_root).resolve() if output_root else PROJECT_ROOT
    ROOT = OUTPUT_ROOT
    raw_root = Path(input_root).resolve() if input_root else OUTPUT_ROOT / "raw"
    MANIFEST = Path(manifest).resolve() if manifest else PROJECT_ROOT / "config" / "species_manifest.tsv"
    DATASETS = resolve_executable(datasets or "datasets")
    RAW_BULK = raw_root / "ncbi_bulk"
    RAW_ASSEMBLIES = raw_root / "assembly_packages"
    INDEXES = OUTPUT_ROOT / "indexes"
    ORTHOLOGY = OUTPUT_ROOT / "orthology"
    SELECTION = OUTPUT_ROOT / "selection"
    QC = OUTPUT_ROOT / "qc"
    FASTAS = OUTPUT_ROOT / "fastas"
    REPORTS = OUTPUT_ROOT / "reports"
    LOGS = OUTPUT_ROOT / "logs"
    RESULTS = OUTPUT_ROOT / "results"


def rel_path(path: Path) -> str:
    for base in [OUTPUT_ROOT, PROJECT_ROOT]:
        try:
            return str(path.relative_to(base))
        except ValueError:
            continue
    return str(path)


def resolve_executable(name: str) -> Path:
    path = Path(name)
    if path.exists():
        return path.resolve()
    if path.parent != Path("."):
        return path
    for candidate in [Path.cwd() / "bin" / name, PROJECT_ROOT / "bin" / name]:
        if candidate.exists():
            return candidate.resolve()
    found = shutil.which(name)
    if found:
        return Path(found)
    return path


def ensure_dirs() -> None:
    for path in [
        RAW_BULK,
        RAW_ASSEMBLIES,
        INDEXES,
        ORTHOLOGY,
        SELECTION,
        QC,
        FASTAS,
        REPORTS,
        LOGS,
        RESULTS,
    ]:
        path.mkdir(parents=True, exist_ok=True)


def read_manifest(path: Optional[Path] = None) -> List[Species]:
    if path is None:
        path = MANIFEST
    rows: List[Species] = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {
            "token",
            "taxid",
            "scientific_name",
            "common_name",
            "gcf_accession",
            "outgroup",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest missing columns: {sorted(missing)}")
        for row in reader:
            rows.append(
                Species(
                    token=row["token"],
                    taxid=int(row["taxid"]),
                    scientific_name=row["scientific_name"],
                    common_name=row["common_name"],
                    gcf_accession=row["gcf_accession"],
                    outgroup=row["outgroup"].lower() == "true",
                )
            )
    tokens = [s.token for s in rows]
    taxids = [s.taxid for s in rows]
    if len(rows) < 2:
        raise ValueError(f"Expected at least 2 species, found {len(rows)}")
    if len(tokens) != len(set(tokens)):
        raise ValueError("Manifest token values are not unique")
    if len(taxids) != len(set(taxids)):
        raise ValueError("Manifest taxid values are not unique")
    bad_tokens = [t for t in tokens if re.search(r"\s", t)]
    if bad_tokens:
        raise ValueError(f"Manifest tokens contain whitespace: {bad_tokens}")
    return rows


def run(cmd: Sequence[str], *, retries: int = 1, cwd: Optional[Path] = None) -> None:
    if cwd is None:
        cwd = PROJECT_ROOT
    cmd_text = " ".join(str(x) for x in cmd)
    for attempt in range(1, retries + 1):
        log(f"RUN attempt {attempt}/{retries}: {cmd_text}")
        result = subprocess.run(cmd, cwd=cwd)
        if result.returncode == 0:
            return
        if attempt < retries:
            sleep_for = min(60, 10 * attempt)
            log(f"Command failed with {result.returncode}; sleeping {sleep_for}s")
            time.sleep(sleep_for)
    raise RuntimeError(f"Command failed after {retries} attempts: {cmd_text}")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)
        fh.write("\n")


def maybe_to_parquet_and_tsv(df: pd.DataFrame, parquet: Path, tsv: Optional[Path] = None) -> None:
    parquet.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(parquet, index=False)
    if tsv is not None:
        df.to_csv(tsv, sep="\t", index=False)


def latest_file(pattern: str) -> Optional[Path]:
    files = sorted(PROJECT_ROOT.glob(pattern))
    return files[-1] if files else None


def latest_bulk_file(prefix: str) -> Optional[Path]:
    files = sorted(RAW_BULK.glob(f"{prefix}*"))
    return files[-1] if files else None


def open_text_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("rt", newline="")


def stage0_preflight(species: List[Species], force: bool = False, offline: bool = False) -> None:
    out = REPORTS / "preflight.json"
    if out.exists() and not force:
        log("Stage 0 preflight already exists; skipping")
        return
    if offline:
        rows = []
        for sp in species:
            package_dir = RAW_ASSEMBLIES / sp.token
            rows.append(
                {
                    "token": sp.token,
                    "taxid": sp.taxid,
                    "gcf_accession": sp.gcf_accession,
                    "package_path": rel_path(package_dir),
                    "package_present": package_done(package_dir),
                    "status": "pass" if package_done(package_dir) else "missing_package",
                }
            )
        df = pd.DataFrame(rows)
        df.to_csv(REPORTS / "gcf_release_ids.tsv", sep="\t", index=False)
        status = "pass" if (df["status"] == "pass").all() else "fail"
        write_json(
            out,
            {
                "status": status,
                "mode": "offline",
                "checked_at": datetime.now().isoformat(),
                "species_count": len(species),
                "tokens": [s.token for s in species],
                "gcf_records": rows,
            },
        )
        if status != "pass":
            raise RuntimeError("Offline preflight failed; see reports/preflight.json")
        log("Stage 0 offline preflight complete")
        return
    if not DATASETS.exists():
        raise FileNotFoundError(f"NCBI datasets CLI not found: {DATASETS}")

    rows: List[Dict[str, Any]] = []
    for sp in species:
        cmd = [
            str(DATASETS),
            "summary",
            "genome",
            "accession",
            sp.gcf_accession,
            "--as-json-lines",
        ]
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"datasets summary failed for {sp.gcf_accession}: {proc.stderr}")
        record = None
        for line in proc.stdout.splitlines():
            if line.strip():
                record = json.loads(line)
                break
        if record is None:
            raise RuntimeError(f"No summary record for {sp.gcf_accession}")
        assembly_info = record.get("assemblyInfo") or record.get("assembly_info") or {}
        annotation_info = record.get("annotationInfo") or record.get("annotation_info") or {}
        rows.append(
            {
                "token": sp.token,
                "taxid": sp.taxid,
                "gcf_accession": sp.gcf_accession,
                "summary_accession": record.get("accession"),
                "assembly_name": assembly_info.get("assemblyName")
                or assembly_info.get("assembly_name"),
                "assembly_level": assembly_info.get("assemblyLevel")
                or assembly_info.get("assembly_level"),
                "annotation_release_id": annotation_info.get("name")
                or annotation_info.get("release_id")
                or annotation_info.get("releaseId"),
                "annotation_release_date": annotation_info.get("releaseDate")
                or annotation_info.get("release_date"),
                "organism": record.get("organism", {}).get("organism_name"),
                "status": "pass" if record.get("accession") == sp.gcf_accession else "accession_mismatch",
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(REPORTS / "gcf_release_ids.tsv", sep="\t", index=False)
    status = "pass" if (df["status"] == "pass").all() else "fail"
    write_json(
        out,
        {
            "status": status,
            "checked_at": datetime.now().isoformat(),
            "species_count": len(species),
            "tokens": [s.token for s in species],
            "gcf_records": rows,
        },
    )
    if status != "pass":
        raise RuntimeError("Preflight failed; see reports/preflight.json")
    log("Stage 0 preflight complete")


def stage1_bulk(force: bool = False, offline: bool = False) -> Tuple[Path, Path]:
    existing_orth = latest_bulk_file("gene_orthologs")
    existing_info = latest_bulk_file("gene_info")
    if offline:
        if not existing_orth or not existing_info:
            raise FileNotFoundError(
                f"Offline mode requires gene_orthologs* and gene_info* under {RAW_BULK}"
            )
        log("Stage 1 using offline NCBI Gene bulk fixtures")
        return existing_orth, existing_info
    if existing_orth and existing_info and not force:
        log("Stage 1 bulk files already exist; skipping")
        return existing_orth, existing_info

    date = datetime.now().strftime("%Y%m%d")
    orth = RAW_BULK / f"gene_orthologs.{date}.gz"
    info = RAW_BULK / f"gene_info.{date}.gz"
    run(
        [
            "curl",
            "-L",
            "--fail",
            "-o",
            str(orth),
            "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz",
        ],
        retries=3,
    )
    run(
        [
            "curl",
            "-L",
            "--fail",
            "-o",
            str(info),
            "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz",
        ],
        retries=3,
    )
    with (RAW_BULK / "MANIFEST.checksums").open("w") as fh:
        for path in [orth, info]:
            fh.write(f"{sha256_file(path)}  {path.name}\n")
    log("Stage 1 NCBI Gene FTP bulk freeze complete")
    return orth, info


def package_done(package_dir: Path) -> bool:
    def usable(path: Path) -> bool:
        return path.is_file() and not any(part.startswith("._") for part in path.parts)

    return any(usable(p) for p in package_dir.rglob("dataset_catalog.json")) and any(
        usable(p) for p in package_dir.rglob("*.gff")
    )


def stage2_download_assemblies(
    species: List[Species],
    force: bool = False,
    offline: bool = False,
    include: str = "cds,protein,rna,gff3,seq-report",
) -> None:
    if offline:
        missing = [sp.token for sp in species if not package_done(RAW_ASSEMBLIES / sp.token)]
        if missing:
            raise FileNotFoundError(f"Offline mode missing assembly packages: {missing}")
        log("Stage 2 using offline assembly packages")
        return
    checksums: List[Dict[str, Any]] = []
    for sp in species:
        package_dir = RAW_ASSEMBLIES / sp.token
        zip_path = RAW_ASSEMBLIES / f"{sp.token}.zip"
        if package_done(package_dir) and not force:
            log(f"Stage 2 package exists for {sp.token}; skipping")
        else:
            if package_dir.exists() and force:
                shutil.rmtree(package_dir)
            package_dir.mkdir(parents=True, exist_ok=True)
            run(
                [
                    str(DATASETS),
                    "download",
                    "genome",
                    "accession",
                    sp.gcf_accession,
                    "--include",
                    include,
                    "--filename",
                    str(zip_path),
                    "--no-progressbar",
                ],
                retries=3,
            )
            if package_dir.exists():
                shutil.rmtree(package_dir)
            package_dir.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(zip_path) as zf:
                zf.extractall(package_dir)
            log(f"Stage 2 downloaded/extracted {sp.token}")
        if zip_path.exists():
            checksums.append(
                {
                    "token": sp.token,
                    "gcf_accession": sp.gcf_accession,
                    "zip_path": rel_path(zip_path),
                    "sha256": sha256_file(zip_path),
                    "size_bytes": zip_path.stat().st_size,
                }
            )
    pd.DataFrame(checksums).to_csv(REPORTS / "assembly_package_checksums.tsv", sep="\t", index=False)
    log("Stage 2 assembly package downloads complete")


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


def first_attr(attrs: Dict[str, str], names: Sequence[str]) -> Optional[str]:
    for name in names:
        value = attrs.get(name)
        if value:
            return value
    return None


def split_attr_values(value: Optional[str]) -> List[str]:
    if not value:
        return []
    return [v for v in re.split(r"[, ]+", value) if v]


def extract_geneid_from_text(text: str) -> Optional[int]:
    match = GENEID_RE.search(text or "")
    return int(match.group(1)) if match else None


def extract_geneid(attrs: Dict[str, str]) -> Optional[int]:
    for key in ["Dbxref", "db_xref", "Ontology_term"]:
        gid = extract_geneid_from_text(attrs.get(key, ""))
        if gid is not None:
            return gid
    return None


def clean_feature_id(value: Optional[str]) -> Optional[str]:
    if not value:
        return None
    return value.strip()


def accession_from_text(text: str, regex: re.Pattern[str]) -> Optional[str]:
    match = regex.search(text or "")
    return match.group(0) if match else None


def transcript_accession(attrs: Dict[str, str]) -> Optional[str]:
    for name in ["transcript_id", "Name", "ID", "Dbxref", "Parent"]:
        acc = accession_from_text(attrs.get(name, ""), TRANSCRIPT_RE)
        if acc:
            return acc
    return None


def protein_accession(attrs: Dict[str, str]) -> Optional[str]:
    for name in ["protein_id", "Name", "ID", "Dbxref"]:
        acc = accession_from_text(attrs.get(name, ""), PROTEIN_RE)
        if acc:
            return acc
    return None


def accession_class(acc: Optional[str]) -> str:
    if not acc:
        return "unknown"
    prefix = acc.split("_", 1)[0]
    return {
        "NM": "validated_mrna",
        "XM": "predicted_mrna",
        "NR": "validated_ncrna",
        "XR": "predicted_ncrna",
    }.get(prefix, "other")


def select_category_from_attrs(attrs: Dict[str, str]) -> str:
    values: List[str] = []
    for key in ["tag", "Note", "exception", "gene_biotype", "gbkey"]:
        raw = attrs.get(key)
        if raw:
            values.extend(v.strip() for v in raw.split(",") if v.strip())
    lowered = {v.lower(): v for v in values}
    for needle in ["mane select", "refseq select"]:
        if needle in lowered:
            return lowered[needle]
    for v in values:
        if "select" in v.lower():
            return v
    return ""


def parse_bracket_fields(description: str) -> Dict[str, str]:
    fields: Dict[str, str] = {}
    for key, value in re.findall(r"\[([^=\]]+)=([^\]]*)\]", description):
        fields[key] = value
    return fields


def fasta_lengths_and_headers(path: Optional[Path]) -> Tuple[Dict[str, int], Dict[str, str], Dict[str, str]]:
    lengths: Dict[str, int] = {}
    headers: Dict[str, str] = {}
    sequences: Dict[str, str] = {}
    if not path or not path.exists():
        return lengths, headers, sequences
    for rec in SeqIO.parse(str(path), "fasta"):
        seq = str(rec.seq).upper()
        identifiers = {rec.id}
        bracket = parse_bracket_fields(rec.description)
        for key in ["protein_id", "transcript_id", "gene"]:
            if bracket.get(key):
                identifiers.add(bracket[key])
        prot = accession_from_text(rec.description, PROTEIN_RE)
        tx = accession_from_text(rec.description, TRANSCRIPT_RE)
        if prot:
            identifiers.add(prot)
        if tx:
            identifiers.add(tx)
        for ident in identifiers:
            lengths[ident] = len(seq)
            headers[ident] = rec.description
            sequences[ident] = seq
    return lengths, headers, sequences


def locate_package_files(package_dir: Path) -> Dict[str, Optional[Path]]:
    def usable_file(path: Path) -> bool:
        return path.is_file() and not any(part.startswith("._") for part in path.parts)

    files = [p for p in package_dir.rglob("*") if usable_file(p)]

    def choose(predicate) -> Optional[Path]:
        candidates = sorted([p for p in files if predicate(p)])
        return candidates[0] if candidates else None

    return {
        "catalog": choose(lambda p: p.name == "dataset_catalog.json"),
        "assembly_report": choose(lambda p: p.name == "assembly_data_report.jsonl"),
        "sequence_report": choose(lambda p: p.name == "sequence_report.jsonl"),
        "gff3": choose(lambda p: p.suffix in {".gff", ".gff3"} and "genomic" in p.name),
        "cds": choose(lambda p: p.suffix in {".fna", ".fa", ".fasta"} and "cds" in p.name.lower()),
        "rna": choose(lambda p: p.suffix in {".fna", ".fa", ".fasta"} and "rna" in p.name.lower()),
        "protein": choose(lambda p: p.suffix in {".faa", ".fa", ".fasta"} and ("protein" in p.name.lower() or p.suffix == ".faa")),
    }


def parse_gene_info(gene_info_path: Path, target_taxids: set[int]) -> Dict[Tuple[int, int], Dict[str, str]]:
    out: Dict[Tuple[int, int], Dict[str, str]] = {}
    with open_text_maybe_gzip(gene_info_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames and reader.fieldnames[0].startswith("#"):
            reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
        for row in reader:
            try:
                taxid = int(row["tax_id"])
                if taxid not in target_taxids:
                    continue
                geneid = int(row["GeneID"])
            except Exception:
                continue
            out[(taxid, geneid)] = {
                "symbol": row.get("Symbol", ""),
                "chromosome": row.get("chromosome", ""),
                "type_of_gene": row.get("type_of_gene", ""),
                "description": row.get("description", ""),
            }
    return out


def read_first_jsonl(path: Optional[Path]) -> Dict[str, Any]:
    if not path or not path.exists():
        return {}
    with path.open() as fh:
        for line in fh:
            if line.strip():
                return json.loads(line)
    return {}


def parse_gff3(
    gff_path: Path,
    token: str,
    taxid: int,
    gene_info: Dict[Tuple[int, int], Dict[str, str]],
    geneid_filter: Optional[set[int]] = None,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    gene_rows: List[Dict[str, Any]] = []
    transcript_rows: List[Dict[str, Any]] = []
    cds_rows_raw: List[Dict[str, Any]] = []
    feature_to_geneid: Dict[str, int] = {}
    feature_to_symbol: Dict[str, str] = {}
    transcript_to_geneid: Dict[str, int] = {}
    transcript_to_acc: Dict[str, str] = {}
    transcript_to_select: Dict[str, str] = {}

    with gff_path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attr_text = parts
            attrs = parse_attrs(attr_text)
            fid = clean_feature_id(attrs.get("ID"))
            gid = extract_geneid(attrs)
            symbol = first_attr(attrs, ["gene", "Name", "gene_synonym"]) or ""

            if ftype == "gene":
                if gid is None:
                    continue
                if geneid_filter is not None and gid not in geneid_filter:
                    continue
                if fid:
                    feature_to_geneid[fid] = gid
                    feature_to_symbol[fid] = symbol
                info = gene_info.get((taxid, gid), {})
                gene_rows.append(
                    {
                        "token": token,
                        "taxid": taxid,
                        "GeneID": gid,
                        "symbol": symbol or info.get("symbol", ""),
                        "chromosome": info.get("chromosome", seqid),
                        "type_of_gene": info.get("type_of_gene", attrs.get("gene_biotype", "")),
                        "gene_biotype": attrs.get("gene_biotype", ""),
                        "feature_id": fid or "",
                        "gff3_seqid": seqid,
                        "start": int(start),
                        "end": int(end),
                        "strand": strand,
                        "source": "assembly_exact",
                    }
                )
                continue

            parent_ids = split_attr_values(attrs.get("Parent"))
            parent_geneid = gid
            for parent in parent_ids:
                if parent_geneid is None and parent in feature_to_geneid:
                    parent_geneid = feature_to_geneid[parent]
                if not symbol and parent in feature_to_symbol:
                    symbol = feature_to_symbol[parent]

            if ftype in {"mRNA", "transcript", "primary_transcript"} or ftype.endswith("RNA"):
                if parent_geneid is None:
                    continue
                if geneid_filter is not None and parent_geneid not in geneid_filter:
                    continue
                tx_acc = transcript_accession(attrs)
                select_category = select_category_from_attrs(attrs)
                if fid:
                    transcript_to_geneid[fid] = parent_geneid
                    if tx_acc:
                        transcript_to_acc[fid] = tx_acc
                    transcript_to_select[fid] = select_category
                transcript_rows.append(
                    {
                        "token": token,
                        "taxid": taxid,
                        "GeneID": parent_geneid,
                        "transcript_accession": tx_acc or "",
                        "transcript_feature_id": fid or "",
                        "transcript_type": ftype,
                        "select_category": select_category,
                        "transcript_length": int(end) - int(start) + 1,
                        "gff3_seqid": seqid,
                        "start": int(start),
                        "end": int(end),
                        "strand": strand,
                        "source": "assembly_exact",
                    }
                )
                continue

            if ftype == "CDS":
                prot = protein_accession(attrs)
                partial = any(
                    [
                        "partial=true" in attr_text.lower(),
                        attrs.get("partial", "").lower() == "true",
                        bool(attrs.get("start_range")),
                        bool(attrs.get("end_range")),
                    ]
                )
                for parent in parent_ids or [""]:
                    row_gid = parent_geneid
                    if row_gid is None and parent in transcript_to_geneid:
                        row_gid = transcript_to_geneid[parent]
                    if row_gid is None:
                        continue
                    if geneid_filter is not None and row_gid not in geneid_filter:
                        continue
                    tx_acc = transcript_to_acc.get(parent) or transcript_accession(attrs) or ""
                    cds_rows_raw.append(
                        {
                            "token": token,
                            "taxid": taxid,
                            "GeneID": row_gid,
                            "transcript_accession": tx_acc,
                            "transcript_feature_id": parent,
                            "cds_accession": prot or f"CDS:{seqid}:{start}-{end}:{strand}:{parent}",
                            "protein_accession": prot or "",
                            "select_category": transcript_to_select.get(parent, ""),
                            "is_partial": partial,
                            "source": "assembly_exact",
                        }
                    )

    # CDS features are exon-level in GFF3; deduplicate to one coding product row.
    dedup: Dict[Tuple[int, str, str, str], Dict[str, Any]] = {}
    for row in cds_rows_raw:
        key = (
            int(row["GeneID"]),
            str(row["transcript_accession"]),
            str(row["cds_accession"]),
            str(row["protein_accession"]),
        )
        if key not in dedup:
            dedup[key] = row
        else:
            dedup[key]["is_partial"] = bool(dedup[key]["is_partial"] or row["is_partial"])
            if not dedup[key].get("select_category") and row.get("select_category"):
                dedup[key]["select_category"] = row["select_category"]
    return gene_rows, transcript_rows, list(dedup.values())


def build_indexes_for_species(
    sp: Species,
    gene_info: Dict[Tuple[int, int], Dict[str, str]],
    geneid_filter: Optional[set[int]] = None,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    package_dir = RAW_ASSEMBLIES / sp.token
    files = locate_package_files(package_dir)
    if not files["gff3"] or not files["cds"]:
        raise FileNotFoundError(f"Missing GFF3 or CDS FASTA for {sp.token}: {files}")

    assembly_record = read_first_jsonl(files.get("assembly_report"))
    assembly_info = assembly_record.get("assemblyInfo") or assembly_record.get("assembly_info") or {}
    annotation_info = assembly_record.get("annotationInfo") or assembly_record.get("annotation_info") or {}
    assembly_rows = [
        {
            "token": sp.token,
            "taxid": sp.taxid,
            "gcf_accession": sp.gcf_accession,
            "summary_accession": assembly_record.get("accession", ""),
            "assembly_name": assembly_info.get("assemblyName")
            or assembly_info.get("assembly_name", ""),
            "annotation_release_id": annotation_info.get("name")
            or annotation_info.get("release_id")
            or annotation_info.get("releaseId", ""),
            "annotation_release_date": annotation_info.get("releaseDate")
            or annotation_info.get("release_date", ""),
            "package_path": rel_path(package_dir),
            "gff3_path": rel_path(files["gff3"]) if files["gff3"] else "",
            "cds_fasta_path": rel_path(files["cds"]) if files["cds"] else "",
            "rna_fasta_path": rel_path(files["rna"]) if files["rna"] else "",
            "protein_fasta_path": rel_path(files["protein"]) if files["protein"] else "",
            "seq_report_path": rel_path(files["sequence_report"]) if files["sequence_report"] else "",
        }
    ]

    rna_lengths, rna_headers, _rna_sequences = fasta_lengths_and_headers(files.get("rna"))
    protein_lengths, protein_headers, _protein_sequences = fasta_lengths_and_headers(files.get("protein"))
    _cds_lengths, cds_headers, cds_sequences = fasta_lengths_and_headers(files.get("cds"))

    gene_rows, transcript_rows, cds_rows = parse_gff3(
        files["gff3"], sp.token, sp.taxid, gene_info, geneid_filter
    )

    # Fill transcript spliced lengths from RNA FASTA when available.
    for row in transcript_rows:
        acc = row.get("transcript_accession")
        if acc and acc in rna_lengths:
            row["transcript_length"] = rna_lengths[acc]

    # Build CDS rows by joining GFF3 product rows to CDS FASTA by protein accession.
    final_cds_rows: List[Dict[str, Any]] = []
    seen_cds_keys: set[Tuple[int, str, str]] = set()
    for row in cds_rows:
        protein_acc = row.get("protein_accession") or ""
        sequence = cds_sequences.get(protein_acc, "")
        header = cds_headers.get(protein_acc, "")
        if not sequence:
            # Some CDS FASTA IDs contain the protein accession in the description
            for ident, seq in cds_sequences.items():
                hdr = cds_headers.get(ident, "")
                if protein_acc and protein_acc in hdr:
                    sequence = seq
                    header = hdr
                    break
        cds_len = len(sequence) if sequence else 0
        protein_len = protein_lengths.get(protein_acc, 0)
        if not protein_len and protein_acc:
            for ident, length in protein_lengths.items():
                if ident == protein_acc or protein_acc in protein_headers.get(ident, ""):
                    protein_len = length
                    break
        key = (int(row["GeneID"]), str(row.get("transcript_accession", "")), str(row.get("protein_accession", "")))
        if key in seen_cds_keys:
            continue
        seen_cds_keys.add(key)
        final_cds_rows.append(
            {
                **row,
                "cds_length": cds_len,
                "protein_length": protein_len,
                "cds_fasta_header": header,
                "cds_sequence": sequence,
                "accession_class": accession_class(row.get("transcript_accession", "")),
            }
        )

    # Header-only fallback for CDS FASTA records not connected through GFF3.
    connected_proteins = {r.get("protein_accession") for r in final_cds_rows if r.get("protein_accession")}
    for ident, seq in cds_sequences.items():
        header = cds_headers.get(ident, "")
        bracket = parse_bracket_fields(header)
        protein_acc = bracket.get("protein_id") or accession_from_text(header, PROTEIN_RE) or ident
        if protein_acc in connected_proteins:
            continue
        gid = extract_geneid_from_text(header)
        if gid is None:
            continue
        if geneid_filter is not None and gid not in geneid_filter:
            continue
        tx_acc = bracket.get("transcript_id") or accession_from_text(header, TRANSCRIPT_RE) or ""
        final_cds_rows.append(
            {
                "token": sp.token,
                "taxid": sp.taxid,
                "GeneID": gid,
                "transcript_accession": tx_acc,
                "transcript_feature_id": "",
                "cds_accession": protein_acc,
                "protein_accession": protein_acc if PROTEIN_RE.match(protein_acc) else "",
                "select_category": "",
                "is_partial": "partial" in header.lower(),
                "source": "assembly_exact_header_only",
                "cds_length": len(seq),
                "protein_length": protein_lengths.get(protein_acc, 0),
                "cds_fasta_header": header,
                "cds_sequence": seq,
                "accession_class": accession_class(tx_acc),
            }
        )

    protein_rows = []
    for ident, length in protein_lengths.items():
        if not accession_from_text(ident, PROTEIN_RE):
            continue
        protein_rows.append(
            {
                "token": sp.token,
                "taxid": sp.taxid,
                "protein_accession": ident,
                "protein_length": length,
                "protein_fasta_header": protein_headers.get(ident, ""),
                "source": "assembly_exact",
            }
        )

    return assembly_rows, gene_rows, transcript_rows, final_cds_rows, protein_rows


def stage3_build_indexes(species: List[Species], gene_info_path: Path, force: bool = False) -> None:
    outputs = [
        INDEXES / "assembly_index.parquet",
        INDEXES / "gene_index.parquet",
        INDEXES / "transcript_index.parquet",
        INDEXES / "cds_index.parquet",
        INDEXES / "protein_index.parquet",
    ]
    if all(p.exists() for p in outputs) and not force:
        log("Stage 3 indexes already exist; skipping")
        return
    target_taxids = {s.taxid for s in species}
    log("Parsing gene_info.gz for target taxids")
    gene_info = parse_gene_info(gene_info_path, target_taxids)
    all_assembly: List[Dict[str, Any]] = []
    all_gene: List[Dict[str, Any]] = []
    all_tx: List[Dict[str, Any]] = []
    all_cds: List[Dict[str, Any]] = []
    all_protein: List[Dict[str, Any]] = []
    for sp in species:
        log(f"Indexing assembly package: {sp.token}")
        geneid_filter = (
            GENEID_FILTER_BY_TAXID.get(sp.taxid)
            if GENEID_FILTER_BY_TAXID is not None
            else None
        )
        assembly, gene, tx, cds, protein = build_indexes_for_species(sp, gene_info, geneid_filter)
        all_assembly.extend(assembly)
        all_gene.extend(gene)
        all_tx.extend(tx)
        all_cds.extend(cds)
        all_protein.extend(protein)
        log(
            f"Indexed {sp.token}: genes={len(gene):,} transcripts={len(tx):,} cds={len(cds):,}"
        )

    assembly_df = pd.DataFrame(all_assembly)
    gene_df = pd.DataFrame(all_gene).drop_duplicates(["taxid", "GeneID"], keep="first")
    tx_df = pd.DataFrame(all_tx).drop_duplicates(
        ["taxid", "GeneID", "transcript_accession", "transcript_feature_id"], keep="first"
    )
    cds_df = pd.DataFrame(all_cds)
    if not cds_df.empty:
        cds_df = cds_df.drop_duplicates(
            ["taxid", "GeneID", "transcript_accession", "protein_accession", "cds_accession"],
            keep="first",
        )
    protein_df = pd.DataFrame(all_protein).drop_duplicates(["taxid", "protein_accession"], keep="first")

    maybe_to_parquet_and_tsv(assembly_df, INDEXES / "assembly_index.parquet", INDEXES / "assembly_index.tsv")
    maybe_to_parquet_and_tsv(gene_df, INDEXES / "gene_index.parquet", INDEXES / "gene_index.tsv")
    maybe_to_parquet_and_tsv(tx_df, INDEXES / "transcript_index.parquet", INDEXES / "transcript_index.tsv")
    maybe_to_parquet_and_tsv(cds_df, INDEXES / "cds_index.parquet", INDEXES / "cds_index.tsv")
    maybe_to_parquet_and_tsv(protein_df, INDEXES / "protein_index.parquet", INDEXES / "protein_index.tsv")
    log("Stage 3 assembly-exact indexes complete")


def stage4_orthology_edges(
    species: List[Species],
    orth_path: Path,
    force: bool = False,
    reference_taxid_filter: Optional[int] = None,
) -> None:
    out = ORTHOLOGY / "ortholog_edges.parquet"
    if out.exists() and not force:
        log("Stage 4 orthology edges already exist; skipping")
        return
    target_taxids = {s.taxid for s in species}
    rows: List[Dict[str, Any]] = []
    if reference_taxid_filter is not None:
        log(
            "Parsing gene_orthologs.gz for edges between "
            f"reference taxid {reference_taxid_filter} and {len(target_taxids)} locked taxa"
        )
    else:
        log(f"Parsing gene_orthologs.gz for edges internal to {len(target_taxids)} locked taxa")
    with open_text_maybe_gzip(orth_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames and reader.fieldnames[0].startswith("#"):
            reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
        for row in reader:
            if row.get("relationship") != "Ortholog":
                continue
            try:
                tax_a = int(row["tax_id"])
                gene_a = int(row["GeneID"])
                tax_b = int(row["Other_tax_id"])
                gene_b = int(row["Other_GeneID"])
            except Exception:
                continue
            if reference_taxid_filter is not None:
                if not (
                    (tax_a == reference_taxid_filter and tax_b in target_taxids)
                    or (tax_b == reference_taxid_filter and tax_a in target_taxids)
                ):
                    continue
            elif tax_a not in target_taxids or tax_b not in target_taxids:
                continue
            left = (tax_a, gene_a)
            right = (tax_b, gene_b)
            if right < left:
                left, right = right, left
            rows.append(
                {
                    "taxid_a": left[0],
                    "GeneID_a": left[1],
                    "taxid_b": right[0],
                    "GeneID_b": right[1],
                    "relationship": "Ortholog",
                    "source_snapshot": orth_path.name,
                }
            )
    df = pd.DataFrame(rows).drop_duplicates()
    maybe_to_parquet_and_tsv(df, out, ORTHOLOGY / "ortholog_edges.tsv")
    log(f"Stage 4 orthology edges complete: {len(df):,} edges")


def node_name(taxid: int, geneid: int) -> str:
    return f"{taxid}:{geneid}"


def parse_node(node: str) -> Tuple[int, int]:
    taxid, geneid = node.split(":", 1)
    return int(taxid), int(geneid)


def existing_orthology_mode(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    try:
        df = pd.read_parquet(path, columns=["orthology_mode"])
    except Exception:
        return None
    if df.empty or "orthology_mode" not in df.columns:
        return None
    modes = {str(x) for x in df["orthology_mode"].dropna().unique()}
    if len(modes) == 1:
        return next(iter(modes))
    return "mixed"


def stage5_strict_singletons(
    species: List[Species],
    reference_taxid: int,
    force: bool = False,
) -> None:
    passed_out = ORTHOLOGY / "strict_singleton.parquet"
    rejected_out = ORTHOLOGY / "rejected.parquet"
    if passed_out.exists() and rejected_out.exists() and not force:
        if existing_orthology_mode(passed_out) == "strict_singleton":
            log("Stage 5 strict singleton outputs already exist; skipping")
            return
        log("Stage 5 existing orthology outputs use a different mode; rebuilding strict singleton outputs")
    edges = pd.read_parquet(ORTHOLOGY / "ortholog_edges.parquet")
    gene_index = pd.read_parquet(INDEXES / "gene_index.parquet")
    token_by_taxid = {s.taxid: s.token for s in species}
    all_taxids = {s.taxid for s in species}
    expected_species_count = len(species)
    if reference_taxid not in all_taxids:
        raise ValueError(
            f"Reference taxid {reference_taxid} is not present in config/species_manifest.tsv"
        )

    gene_meta: Dict[Tuple[int, int], Dict[str, Any]] = {}
    for row in gene_index.to_dict("records"):
        gene_meta[(int(row["taxid"]), int(row["GeneID"]))] = row

    graph = nx.Graph()
    for row in edges.itertuples(index=False):
        a = node_name(int(row.taxid_a), int(row.GeneID_a))
        b = node_name(int(row.taxid_b), int(row.GeneID_b))
        graph.add_edge(a, b)

    passed_rows: List[Dict[str, Any]] = []
    rejected_rows: List[Dict[str, Any]] = []
    component_rows: List[Dict[str, Any]] = []
    comp_idx = 0
    for component in nx.connected_components(graph):
        comp_idx += 1
        nodes = sorted(component)
        parsed = [parse_node(n) for n in nodes]
        tax_counts = Counter(t for t, _g in parsed)
        gene_count = len(parsed)
        reference_genes = [g for t, g in parsed if t == reference_taxid]
        component_rows.append(
            {
                "component_id": comp_idx,
                "gene_count": gene_count,
                "taxid_count": len(tax_counts),
                "has_reference": bool(reference_genes),
                "nodes": ";".join(nodes),
            }
        )

        reason = ""
        if not reference_genes:
            reason = "missing_reference"
        elif len(reference_genes) != 1:
            reason = "multiple_reference_genes"
        elif gene_count != expected_species_count:
            reason = "component_size_not_expected_species_count"
        elif set(tax_counts) != all_taxids:
            reason = "missing_taxon"
        elif any(v != 1 for v in tax_counts.values()):
            reason = "multi_gene_per_taxon"
        else:
            for taxid, geneid in parsed:
                meta = gene_meta.get((taxid, geneid))
                if meta is None:
                    reason = "missing_assembly_gene"
                    break
                chrom = str(meta.get("chromosome", ""))
                gene_type = str(meta.get("type_of_gene", ""))
                if chrom.upper() in {"MT", "MITOCHONDRION"}:
                    reason = "mt_gene"
                    break
                if gene_type and gene_type != "protein-coding":
                    reason = "non_protein_coding"
                    break

        if reason:
            rejected_rows.append(
                {
                    "component_id": comp_idx,
                    "reason": reason,
                    "gene_count": gene_count,
                    "taxid_count": len(tax_counts),
                    "reference_taxid": reference_taxid,
                    "reference_GeneIDs": ",".join(str(g) for g in reference_genes),
                    "nodes": ";".join(nodes),
                }
            )
            continue

        reference_geneid = reference_genes[0]
        reference_meta = gene_meta[(reference_taxid, reference_geneid)]
        reference_symbol = reference_meta.get("symbol") or f"GeneID_{reference_geneid}"
        family_id = f"{reference_symbol}__{reference_geneid}"
        for taxid, geneid in sorted(parsed):
            meta = gene_meta[(taxid, geneid)]
            passed_rows.append(
                {
                    "family_id": family_id,
                    "component_id": comp_idx,
                    "orthology_mode": "strict_singleton",
                    "family_expected_count": expected_species_count,
                    "family_min_sequences": expected_species_count,
                    "query_symbol": reference_symbol,
                    "reference_taxid": reference_taxid,
                    "reference_token": token_by_taxid[reference_taxid],
                    "reference_GeneID": reference_geneid,
                    "reference_symbol": reference_symbol,
                    # Backward-compatible aliases. These are reference fields,
                    # not necessarily human fields when --reference-taxid != 9606.
                    "human_GeneID": reference_geneid,
                    "human_symbol": reference_symbol,
                    "taxid": taxid,
                    "token": token_by_taxid[taxid],
                    "GeneID": geneid,
                    "species_symbol": meta.get("symbol", ""),
                }
            )

    comp_df = pd.DataFrame(component_rows)
    passed_df = pd.DataFrame(passed_rows)
    rejected_df = pd.DataFrame(rejected_rows)
    maybe_to_parquet_and_tsv(comp_df, ORTHOLOGY / "candidate_components.parquet", ORTHOLOGY / "candidate_components.tsv")
    maybe_to_parquet_and_tsv(passed_df, passed_out, ORTHOLOGY / "strict_singleton.tsv")
    maybe_to_parquet_and_tsv(rejected_df, rejected_out, ORTHOLOGY / "rejected.tsv")
    log(
        f"Stage 5 strict singletons complete: families={passed_df['family_id'].nunique() if not passed_df.empty else 0:,}; "
        f"rejected_components={len(rejected_df):,}"
    )


def reference_query_dir(symbol: str, geneid: int) -> Path:
    base = re.sub(r"[^A-Za-z0-9_.-]+", "_", (symbol or f"GeneID_{geneid}").strip()).strip("._")
    return RESULTS / "reference_gene_1to1" / (base or f"GeneID_{geneid}")


def validate_gene_meta(meta: Optional[Dict[str, Any]], *, role: str) -> str:
    if meta is None:
        return f"{role}_gene_not_in_gcf_annotation"
    chrom = str(meta.get("chromosome", ""))
    gene_type = str(meta.get("type_of_gene", ""))
    if chrom.upper() in {"MT", "MITOCHONDRION"}:
        return "mitochondrial_gene"
    if gene_type and gene_type != "protein-coding":
        return "non_protein_coding"
    return ""


def resolve_reference_gene(
    gene_index: pd.DataFrame,
    reference_taxid: int,
    reference_symbol: Optional[str],
    reference_gene_id: Optional[int],
) -> Dict[str, Any]:
    ref_genes = gene_index[gene_index["taxid"].astype(int) == int(reference_taxid)].copy()
    if ref_genes.empty:
        raise RuntimeError(f"No genes indexed for reference taxid {reference_taxid}")

    if reference_gene_id is not None:
        matches = ref_genes[ref_genes["GeneID"].astype(int) == int(reference_gene_id)]
        method = "geneid"
        query_symbol = reference_symbol or ""
    elif reference_symbol:
        query_symbol = reference_symbol
        matches = ref_genes[ref_genes["symbol"].fillna("") == reference_symbol]
        method = "symbol_exact"
        if matches.empty:
            ci = ref_genes[ref_genes["symbol"].fillna("").str.lower() == reference_symbol.lower()]
            if len(ci) == 1:
                matches = ci
                method = "symbol_case_insensitive_unique"
    else:
        raise ValueError("reference_gene_1to1_present_species mode requires --reference-symbol or --reference-gene-id")

    if matches.empty:
        raise RuntimeError(f"Reference gene not found for taxid {reference_taxid}: {reference_symbol or reference_gene_id}")
    if len(matches) > 1:
        genes = ",".join(str(int(x)) for x in matches["GeneID"].tolist())
        raise RuntimeError(
            f"Reference symbol is ambiguous in taxid {reference_taxid}: {reference_symbol}; GeneIDs={genes}. "
            "Use --reference-gene-id."
        )

    row = matches.iloc[0].to_dict()
    reason = validate_gene_meta(row, role="reference")
    if reason:
        raise RuntimeError(f"Reference gene failed validation: {reason}")
    return {
        "reference_taxid": int(reference_taxid),
        "query_symbol": query_symbol,
        "resolved_symbol": row.get("symbol") or f"GeneID_{int(row['GeneID'])}",
        "reference_gene_id": int(row["GeneID"]),
        "gene_type": row.get("type_of_gene", ""),
        "chromosome": row.get("chromosome", ""),
        "resolution_method": method,
        "status": "pass",
    }


def edge_neighbors(
    edges: pd.DataFrame,
    taxid: int,
    geneid: int,
    other_taxid: Optional[int] = None,
) -> set[int]:
    out: set[int] = set()
    left = edges[(edges["taxid_a"].astype(int) == int(taxid)) & (edges["GeneID_a"].astype(int) == int(geneid))]
    if other_taxid is not None:
        left = left[left["taxid_b"].astype(int) == int(other_taxid)]
    out.update(int(x) for x in left["GeneID_b"].tolist())
    right = edges[(edges["taxid_b"].astype(int) == int(taxid)) & (edges["GeneID_b"].astype(int) == int(geneid))]
    if other_taxid is not None:
        right = right[right["taxid_a"].astype(int) == int(other_taxid)]
    out.update(int(x) for x in right["GeneID_a"].tolist())
    return out


def stage5_reference_gene_1to1_present_species(
    species: List[Species],
    reference_taxid: int,
    reference_symbol: Optional[str],
    reference_gene_id: Optional[int],
    include_reference: bool,
    min_sequences: int,
    force: bool = False,
) -> None:
    passed_out = ORTHOLOGY / "strict_singleton.parquet"
    rejected_out = ORTHOLOGY / "rejected.parquet"
    if passed_out.exists() and rejected_out.exists() and not force:
        if existing_orthology_mode(passed_out) == "reference_gene_1to1_present_species":
            log("Stage 5 reference-gene 1:1 outputs already exist; skipping")
            return
        log("Stage 5 existing orthology outputs use a different mode; rebuilding reference-gene 1:1 outputs")

    edges = pd.read_parquet(ORTHOLOGY / "ortholog_edges.parquet")
    gene_index = pd.read_parquet(INDEXES / "gene_index.parquet")
    token_by_taxid = {s.taxid: s.token for s in species}
    species_by_taxid = {s.taxid: s for s in species}
    all_taxids = {s.taxid for s in species}
    if reference_taxid not in all_taxids:
        raise ValueError(f"Reference taxid {reference_taxid} is not present in the manifest")

    gene_meta: Dict[Tuple[int, int], Dict[str, Any]] = {}
    for row in gene_index.to_dict("records"):
        gene_meta[(int(row["taxid"]), int(row["GeneID"]))] = row

    ref = resolve_reference_gene(gene_index, reference_taxid, reference_symbol, reference_gene_id)
    reference_geneid = int(ref["reference_gene_id"])
    reference_symbol_resolved = str(ref["resolved_symbol"])
    family_id = f"{reference_symbol_resolved}__{reference_geneid}__reference_1to1"
    query_dir = reference_query_dir(reference_symbol_resolved, reference_geneid)
    query_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame([ref]).to_csv(query_dir / "reference_gene.tsv", sep="\t", index=False)

    candidate_rows: List[Dict[str, Any]] = []
    rejected_rows: List[Dict[str, Any]] = []
    pass_graph_rows: List[Dict[str, Any]] = []
    passed_rows: List[Dict[str, Any]] = []

    if include_reference:
        ref_meta = gene_meta.get((reference_taxid, reference_geneid))
        reason = validate_gene_meta(ref_meta, role="reference")
        if reason:
            raise RuntimeError(f"Reference gene failed annotation validation: {reason}")
        passed_rows.append(
            {
                "family_id": family_id,
                "component_id": 1,
                "orthology_mode": "reference_gene_1to1_present_species",
                "family_expected_count": 0,  # filled below after target classification
                "family_min_sequences": max(1, int(min_sequences)),
                "reference_included": bool(include_reference),
                "query_symbol": ref["query_symbol"] or reference_symbol_resolved,
                "reference_taxid": reference_taxid,
                "reference_token": token_by_taxid[reference_taxid],
                "reference_GeneID": reference_geneid,
                "reference_symbol": reference_symbol_resolved,
                "human_GeneID": reference_geneid,
                "human_symbol": reference_symbol_resolved,
                "taxid": reference_taxid,
                "token": token_by_taxid[reference_taxid],
                "GeneID": reference_geneid,
                "species_symbol": ref_meta.get("symbol", reference_symbol_resolved),
                "orthology_status": "reference_included",
                "reject_reason": "",
            }
        )

    for sp in species:
        if sp.taxid == reference_taxid:
            continue
        target_genes = edge_neighbors(edges, reference_taxid, reference_geneid, sp.taxid)
        row_base = {
            "reference_symbol": reference_symbol_resolved,
            "reference_gene_id": reference_geneid,
            "target_taxid": sp.taxid,
            "target_species_token": sp.token,
            "target_gene_id": "",
            "target_gene_count_for_reference": len(target_genes),
            "reference_gene_count_for_target": "",
            "edge_direction_summary": "undirected_gene_orthologs",
        }
        if len(target_genes) == 0:
            rejected_rows.append({**row_base, "orthology_status": "reject", "reject_reason": "no_ortholog_edge"})
            candidate_rows.append({**row_base, "orthology_status": "reject", "reject_reason": "no_ortholog_edge"})
            continue

        target_ref_counts: Dict[int, set[int]] = {
            target_gene: edge_neighbors(edges, sp.taxid, target_gene, reference_taxid)
            for target_gene in target_genes
        }
        if len(target_genes) > 1:
            reason = "one_reference_to_many_targets"
            if any(len(refs) > 1 for refs in target_ref_counts.values()):
                reason = "many_to_many_ambiguous"
            rejected_rows.append(
                {
                    **row_base,
                    "target_gene_id": ",".join(str(g) for g in sorted(target_genes)),
                    "reference_gene_count_for_target": ",".join(
                        f"{g}:{len(target_ref_counts[g])}" for g in sorted(target_ref_counts)
                    ),
                    "orthology_status": "reject",
                    "reject_reason": reason,
                }
            )
            candidate_rows.append(rejected_rows[-1])
            continue

        target_geneid = next(iter(target_genes))
        reference_genes_for_target = target_ref_counts[target_geneid]
        row_base["target_gene_id"] = target_geneid
        row_base["reference_gene_count_for_target"] = len(reference_genes_for_target)
        if reference_genes_for_target != {reference_geneid}:
            reason = "many_reference_to_one_target_ambiguous" if len(reference_genes_for_target) > 1 else "inconsistent_ortholog_edge"
            rejected_rows.append({**row_base, "orthology_status": "reject", "reject_reason": reason})
            candidate_rows.append(rejected_rows[-1])
            continue

        meta = gene_meta.get((sp.taxid, target_geneid))
        reason = validate_gene_meta(meta, role="target")
        if reason:
            rejected_rows.append({**row_base, "orthology_status": "reject", "reject_reason": reason})
            candidate_rows.append(rejected_rows[-1])
            continue

        pass_row = {**row_base, "orthology_status": "pass_graph_1to1", "reject_reason": ""}
        candidate_rows.append(pass_row)
        pass_graph_rows.append(pass_row)
        passed_rows.append(
            {
                "family_id": family_id,
                "component_id": 1,
                "orthology_mode": "reference_gene_1to1_present_species",
                "family_expected_count": 0,  # filled below
                "family_min_sequences": max(1, int(min_sequences)),
                "reference_included": bool(include_reference),
                "query_symbol": ref["query_symbol"] or reference_symbol_resolved,
                "reference_taxid": reference_taxid,
                "reference_token": token_by_taxid[reference_taxid],
                "reference_GeneID": reference_geneid,
                "reference_symbol": reference_symbol_resolved,
                "human_GeneID": reference_geneid,
                "human_symbol": reference_symbol_resolved,
                "taxid": sp.taxid,
                "token": token_by_taxid[sp.taxid],
                "GeneID": target_geneid,
                "species_symbol": meta.get("symbol", ""),
                "orthology_status": "pass_graph_1to1",
                "reject_reason": "",
            }
        )

    family_expected_count = len(passed_rows)
    for row in passed_rows:
        row["family_expected_count"] = family_expected_count

    passed_df = pd.DataFrame(passed_rows)
    rejected_df = pd.DataFrame(rejected_rows)
    candidates_df = pd.DataFrame(candidate_rows)
    pass_graph_df = pd.DataFrame(pass_graph_rows)

    component_df = pd.DataFrame(
        [
            {
                "component_id": 1,
                "family_id": family_id,
                "orthology_mode": "reference_gene_1to1_present_species",
                "reference_taxid": reference_taxid,
                "reference_GeneID": reference_geneid,
                "reference_symbol": reference_symbol_resolved,
                "manifest_species_count": len(species),
                "graph_pass_target_species_count": len(pass_graph_df),
                "included_sequence_candidates": len(passed_df),
                "rejected_species_count": len(rejected_df),
            }
        ]
    )

    maybe_to_parquet_and_tsv(component_df, ORTHOLOGY / "candidate_components.parquet", ORTHOLOGY / "candidate_components.tsv")
    maybe_to_parquet_and_tsv(passed_df, passed_out, ORTHOLOGY / "strict_singleton.tsv")
    maybe_to_parquet_and_tsv(rejected_df, rejected_out, ORTHOLOGY / "rejected.tsv")
    candidates_df.to_csv(query_dir / "ortholog_candidates.tsv", sep="\t", index=False)
    pass_graph_df.to_csv(query_dir / "ortholog_pass_graph_1to1.tsv", sep="\t", index=False)
    rejected_df.to_csv(query_dir / "ortholog_rejected_graph.tsv", sep="\t", index=False)
    write_json(
        query_dir / "ortholog_coverage.json",
        {
            "reference_symbol": reference_symbol_resolved,
            "reference_gene_id": reference_geneid,
            "manifest_species_count": len(species),
            "graph_pass_target_species_count": len(pass_graph_df),
            "included_sequence_candidates": len(passed_df),
            "rejected_reasons": rejected_df["reject_reason"].value_counts().to_dict() if not rejected_df.empty else {},
        },
    )
    log(
        "Stage 5 reference-gene 1:1 complete: "
        f"query={reference_symbol_resolved}; included_candidates={len(passed_df):,}; rejected_species={len(rejected_df):,}"
    )


def species_needed_for_reference_gene_mode(
    species: List[Species],
    orth_path: Path,
    info_path: Path,
    reference_taxid: int,
    reference_symbol: Optional[str],
    reference_gene_id: Optional[int],
) -> List[Species]:
    """Return the minimal assembly package set for a reference-gene query.

    For single-gene reference mode, downloading every manifest assembly can be
    unnecessary. NCBI Gene bulk files are enough to identify graph-level 1:1
    target taxa first; assembly packages are then needed only for the reference
    species and those graph-passing target species.
    """
    global GENEID_FILTER_BY_TAXID
    target_taxids = {s.taxid for s in species}
    gene_info = parse_gene_info(info_path, target_taxids)
    gene_df = pd.DataFrame(
        [
            {
                "taxid": taxid,
                "GeneID": geneid,
                **meta,
            }
            for (taxid, geneid), meta in gene_info.items()
        ]
    )
    ref = resolve_reference_gene(gene_df, reference_taxid, reference_symbol, reference_gene_id)
    reference_geneid = int(ref["reference_gene_id"])
    ref_to_targets: Dict[Tuple[int, int], set[int]] = defaultdict(set)
    target_to_refs: Dict[Tuple[int, int], set[int]] = defaultdict(set)
    with open_text_maybe_gzip(orth_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames and reader.fieldnames[0].startswith("#"):
            reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
        for row in reader:
            if row.get("relationship") != "Ortholog":
                continue
            try:
                tax_a = int(row["tax_id"])
                gene_a = int(row["GeneID"])
                tax_b = int(row["Other_tax_id"])
                gene_b = int(row["Other_GeneID"])
            except Exception:
                continue
            if tax_a == reference_taxid and tax_b in target_taxids:
                ref_to_targets[(gene_a, tax_b)].add(gene_b)
                target_to_refs[(tax_b, gene_b)].add(gene_a)
            elif tax_b == reference_taxid and tax_a in target_taxids:
                ref_to_targets[(gene_b, tax_a)].add(gene_a)
                target_to_refs[(tax_a, gene_a)].add(gene_b)

    needed_taxids = {reference_taxid}
    geneid_filter_by_taxid: Dict[int, set[int]] = {reference_taxid: {reference_geneid}}
    for sp in species:
        if sp.taxid == reference_taxid:
            continue
        target_genes = ref_to_targets.get((reference_geneid, sp.taxid), set())
        if len(target_genes) != 1:
            continue
        target_geneid = next(iter(target_genes))
        if target_to_refs.get((sp.taxid, target_geneid), set()) == {reference_geneid}:
            needed_taxids.add(sp.taxid)
            geneid_filter_by_taxid.setdefault(sp.taxid, set()).add(target_geneid)

    scoped = [sp for sp in species if sp.taxid in needed_taxids]
    GENEID_FILTER_BY_TAXID = geneid_filter_by_taxid
    write_json(
        REPORTS / "reference_gene_download_scope.json",
        {
            "reference_symbol": ref["resolved_symbol"],
            "reference_gene_id": reference_geneid,
            "manifest_species_count": len(species),
            "assembly_packages_needed": len(scoped),
            "target_taxids_with_graph_1to1": sorted(needed_taxids - {reference_taxid}),
            "geneids_by_taxid": {
                str(taxid): sorted(geneids)
                for taxid, geneids in sorted(geneid_filter_by_taxid.items())
            },
            "tokens": [sp.token for sp in scoped],
        },
    )
    log(
        "Reference-gene mode assembly download scope: "
        f"{len(scoped):,}/{len(species):,} assemblies "
        f"for {ref['resolved_symbol']} ({reference_geneid})"
    )
    return scoped


def select_rank(select_category: str) -> int:
    text = (select_category or "").lower()
    if "mane select" in text:
        return 0
    if "refseq select" in text:
        return 1
    if "select" in text:
        return 2
    return 3


def accession_rank(cls: str) -> int:
    if cls == "validated_mrna":
        return 0
    if cls == "predicted_mrna":
        return 1
    return 2


def family_completion_status(df: pd.DataFrame, count_col: str = "taxid") -> Dict[str, bool]:
    if df.empty:
        return {}
    out: Dict[str, bool] = {}
    for family_id, fam in df.groupby("family_id", sort=False):
        mode = str(fam["orthology_mode"].iloc[0]) if "orthology_mode" in fam.columns else "strict_singleton"
        expected = int(fam["family_expected_count"].iloc[0]) if "family_expected_count" in fam.columns else fam[count_col].nunique()
        min_sequences = int(fam["family_min_sequences"].iloc[0]) if "family_min_sequences" in fam.columns else expected
        observed = int(fam[count_col].nunique())
        reference_taxid = int(fam["reference_taxid"].iloc[0]) if "reference_taxid" in fam.columns else 9606
        has_reference = ((fam["taxid"].astype(int) == reference_taxid).any()) if "taxid" in fam.columns else True
        reference_required = bool(fam["reference_included"].iloc[0]) if "reference_included" in fam.columns else True
        if mode == "reference_gene_1to1_present_species":
            out[str(family_id)] = bool((has_reference or not reference_required) and observed >= min_sequences)
        else:
            out[str(family_id)] = bool(observed == expected)
    return out


def stage6_select_representatives(expected_species_count: int, force: bool = False) -> None:
    out = SELECTION / "representative_cds.parquet"
    if out.exists() and not force:
        log("Stage 6 representative CDS already exists; skipping")
        return
    singleton = pd.read_parquet(ORTHOLOGY / "strict_singleton.parquet")
    cds = pd.read_parquet(INDEXES / "cds_index.parquet")
    if singleton.empty:
        raise RuntimeError("No ortholog families found; cannot select CDS")
    cds["has_sequence"] = cds["cds_sequence"].fillna("").str.len() > 0
    cds["select_rank"] = cds["select_category"].fillna("").map(select_rank)
    cds["accession_rank"] = cds["accession_class"].fillna("unknown").map(accession_rank)
    cds["protein_length_num"] = pd.to_numeric(cds["protein_length"], errors="coerce").fillna(0).astype(int)
    cds["cds_length_num"] = pd.to_numeric(cds["cds_length"], errors="coerce").fillna(0).astype(int)
    cds["is_partial_bool"] = cds["is_partial"].fillna(False).astype(bool)

    selected: List[Dict[str, Any]] = []
    audit: List[Dict[str, Any]] = []
    grouped = cds[cds["has_sequence"]].groupby(["taxid", "GeneID"], sort=False)
    candidate_map = {key: df.copy() for key, df in grouped}

    for row in singleton.to_dict("records"):
        key = (int(row["taxid"]), int(row["GeneID"]))
        cand = candidate_map.get(key)
        if cand is None or cand.empty:
            audit.append({**row, "status": "fail", "reason": "no_assembly_cds"})
            continue
        usable = cand[~cand["is_partial_bool"]].copy()
        if usable.empty:
            usable = cand.copy()
            partial_fallback = True
        else:
            partial_fallback = False
        usable = usable.sort_values(
            by=[
                "select_rank",
                "protein_length_num",
                "cds_length_num",
                "accession_rank",
                "transcript_accession",
                "protein_accession",
            ],
            ascending=[True, False, False, True, True, True],
            kind="mergesort",
        )
        chosen = usable.iloc[0].to_dict()
        selected.append(
            {
                **row,
                "transcript_accession": chosen.get("transcript_accession", ""),
                "cds_accession": chosen.get("cds_accession", ""),
                "protein_accession": chosen.get("protein_accession", ""),
                "cds_length": int(chosen.get("cds_length_num", 0)),
                "protein_length": int(chosen.get("protein_length_num", 0)),
                "select_category": chosen.get("select_category", ""),
                "accession_class": chosen.get("accession_class", ""),
                "is_partial": bool(chosen.get("is_partial_bool", False)),
                "selection_rule_id": "assembly_exact_select_complete_longest_coding",
                "used_partial_fallback": partial_fallback,
                "cds_sequence": chosen.get("cds_sequence", ""),
                "source": chosen.get("source", "assembly_exact"),
            }
        )
        audit.append(
            {
                **row,
                "status": "pass",
                "reason": "",
                "candidate_count": len(cand),
                "usable_candidate_count": len(usable),
                "used_partial_fallback": partial_fallback,
                "chosen_transcript_accession": chosen.get("transcript_accession", ""),
                "chosen_protein_accession": chosen.get("protein_accession", ""),
            }
        )

    selected_df = pd.DataFrame(selected)
    audit_df = pd.DataFrame(audit)
    if selected_df.empty:
        selected_df = pd.DataFrame(columns=list(singleton.columns) + ["family_selection_complete"])

    # Strict mode requires every locked species. Reference-gene mode is
    # species-level: species without usable CDS are excluded while the query
    # survives if the reference sequence and min sequence count are present.
    complete_map = family_completion_status(selected_df)
    complete_families = {fam for fam, ok in complete_map.items() if ok}
    if "family_id" in selected_df.columns:
        selected_df["family_selection_complete"] = selected_df["family_id"].astype(str).isin(complete_families)
    else:
        selected_df["family_selection_complete"] = False
    maybe_to_parquet_and_tsv(selected_df, out, SELECTION / "representative_cds.tsv")
    maybe_to_parquet_and_tsv(audit_df, SELECTION / "selection_audit.parquet", SELECTION / "selection_audit.tsv")
    log(f"Stage 6 representative CDS complete: complete_families={len(complete_families):,}")


def normalize_cds(seq: str) -> Dict[str, Any]:
    raw = (seq or "").upper().replace("U", "T").replace("-", "")
    result: Dict[str, Any] = {
        "passed": False,
        "fail_reason": "",
        "cds_length_original": len(raw),
        "cds_length_normalized": len(raw),
        "had_terminal_stop": False,
        "terminal_stop_removed": False,
        "n_fraction": 0.0,
        "normalized_sequence": raw,
        "start_codon": raw[:3] if len(raw) >= 3 else "",
    }
    if len(raw) < 60:
        result["fail_reason"] = "length_lt_60"
        return result
    if not DNA_RE.match(raw):
        result["fail_reason"] = "invalid_dna_alphabet"
        return result
    result["n_fraction"] = raw.count("N") / len(raw) if raw else 0.0
    if result["n_fraction"] > 0.01:
        result["fail_reason"] = "n_fraction_gt_0_01"
        return result
    if len(raw) % 3 != 0:
        result["fail_reason"] = "length_not_multiple_of_3"
        return result
    normalized = raw
    terminal = raw[-3:]
    if terminal in STOP_CODONS:
        normalized = raw[:-3]
        result["had_terminal_stop"] = True
        result["terminal_stop_removed"] = True
    for i in range(0, len(normalized), 3):
        codon = normalized[i : i + 3]
        if len(codon) < 3:
            result["fail_reason"] = "normalized_length_not_multiple_of_3"
            return result
        if "N" in codon:
            continue
        if codon in STOP_CODONS:
            result["fail_reason"] = "internal_stop"
            return result
    result["passed"] = True
    result["cds_length_normalized"] = len(normalized)
    result["normalized_sequence"] = normalized
    return result


def stage7_cds_qc(expected_species_count: int, force: bool = False) -> None:
    out = QC / "cds_qc.parquet"
    norm_out = SELECTION / "representative_cds.normalized.parquet"
    if out.exists() and norm_out.exists() and not force:
        log("Stage 7 CDS QC already exists; skipping")
        return
    selected = pd.read_parquet(SELECTION / "representative_cds.parquet")
    selected = selected[selected["family_selection_complete"] == True].copy()
    qc_rows: List[Dict[str, Any]] = []
    norm_rows: List[Dict[str, Any]] = []
    for row in selected.to_dict("records"):
        qc = normalize_cds(row.get("cds_sequence", ""))
        qc_row = {
            "family_id": row["family_id"],
            "human_symbol": row["human_symbol"],
            "taxid": row["taxid"],
            "token": row["token"],
            "GeneID": row["GeneID"],
            "transcript_accession": row.get("transcript_accession", ""),
            "protein_accession": row.get("protein_accession", ""),
            **{k: v for k, v in qc.items() if k != "normalized_sequence"},
        }
        qc_rows.append(qc_row)
        norm_rows.append({**row, **qc})
    qc_df = pd.DataFrame(qc_rows)
    norm_df = pd.DataFrame(norm_rows)
    if not qc_df.empty:
        passed_norm = norm_df[norm_df["passed"] == True].copy()
        complete_map = family_completion_status(passed_norm)
        pass_families = {fam for fam, ok in complete_map.items() if ok}
        qc_df["family_cds_qc_passed"] = qc_df["family_id"].astype(str).isin(pass_families)
        norm_df["family_cds_qc_passed"] = norm_df["family_id"].astype(str).isin(pass_families)
    maybe_to_parquet_and_tsv(qc_df, out, QC / "cds_qc.tsv")
    maybe_to_parquet_and_tsv(norm_df, norm_out, SELECTION / "representative_cds.normalized.tsv")
    pass_family_count = 0
    if not norm_df.empty and "family_cds_qc_passed" in norm_df.columns:
        pass_family_count = norm_df[norm_df["family_cds_qc_passed"]]["family_id"].nunique()
    log(f"Stage 7 CDS QC complete: pass_families={pass_family_count:,}")


def stage8_family_sanity(force: bool = False) -> None:
    out = QC / "family_sanity.parquet"
    if out.exists() and not force:
        log("Stage 8 family sanity already exists; skipping")
        return
    norm = pd.read_parquet(SELECTION / "representative_cds.normalized.parquet")
    if norm.empty or "family_cds_qc_passed" not in norm.columns:
        df = pd.DataFrame(
            columns=[
                "family_id",
                "reference_symbol",
                "species_count",
                "protein_length_min",
                "protein_length_max",
                "protein_length_ratio",
                "cds_length_min",
                "cds_length_max",
                "cds_length_ratio",
                "status",
                "flags",
            ]
        )
        maybe_to_parquet_and_tsv(df, out, QC / "family_sanity.tsv")
        log("Stage 8 family sanity complete: families=0")
        return
    norm = norm[(norm["family_cds_qc_passed"] == True) & (norm["passed"] == True)].copy()
    rows: List[Dict[str, Any]] = []
    for family_id, fam in norm.groupby("family_id"):
        protein_lengths = pd.to_numeric(fam["protein_length"], errors="coerce").fillna(0)
        cds_lengths = pd.to_numeric(fam["cds_length_normalized"], errors="coerce").fillna(0)
        min_prot = int(protein_lengths[protein_lengths > 0].min()) if (protein_lengths > 0).any() else 0
        max_prot = int(protein_lengths.max()) if len(protein_lengths) else 0
        min_cds = int(cds_lengths[cds_lengths > 0].min()) if (cds_lengths > 0).any() else 0
        max_cds = int(cds_lengths.max()) if len(cds_lengths) else 0
        prot_ratio = (max_prot / min_prot) if min_prot else None
        cds_ratio = (max_cds / min_cds) if min_cds else None
        flags = []
        if prot_ratio is not None and prot_ratio > 1.5:
            flags.append("protein_length_ratio_gt_1_5")
        if prot_ratio is not None and prot_ratio > 2.0:
            flags.append("protein_length_ratio_gt_2_0")
        if cds_ratio is not None and cds_ratio > 2.0:
            flags.append("cds_length_ratio_gt_2_0")
        rows.append(
            {
                "family_id": family_id,
                "reference_symbol": fam["reference_symbol"].iloc[0]
                if "reference_symbol" in fam.columns
                else fam["human_symbol"].iloc[0],
                "species_count": fam["taxid"].nunique(),
                "protein_length_min": min_prot,
                "protein_length_max": max_prot,
                "protein_length_ratio": prot_ratio,
                "cds_length_min": min_cds,
                "cds_length_max": max_cds,
                "cds_length_ratio": cds_ratio,
                "status": "flag" if flags else "pass",
                "flags": ",".join(flags),
            }
        )
    df = pd.DataFrame(rows)
    maybe_to_parquet_and_tsv(df, out, QC / "family_sanity.tsv")
    log(f"Stage 8 family sanity complete: families={len(df):,}")


def safe_symbol(symbol: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", symbol.strip())
    cleaned = cleaned.strip("._")
    return cleaned or "UNKNOWN_SYMBOL"


def wrap_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def stage9_write_fastas(species: List[Species], force: bool = False) -> None:
    manifest_out = FASTAS / "manifest.tsv"
    if manifest_out.exists() and not force:
        log("Stage 9 FASTA manifest already exists; skipping")
        return
    norm = pd.read_parquet(SELECTION / "representative_cds.normalized.parquet")
    sanity = pd.read_parquet(QC / "family_sanity.parquet")
    assembly = pd.read_parquet(INDEXES / "assembly_index.parquet")
    assembly_meta = assembly.set_index("token").to_dict("index")
    manifest_columns = [
        "family_id",
        "orthology_mode",
        "reference_symbol",
        "reference_GeneID",
        "reference_taxid",
        "reference_token",
        "human_symbol",
        "human_GeneID",
        "fasta_path",
        "meta_path",
        "sequence_count",
        "sanity_status",
    ]
    keep_families = set(sanity["family_id"]) if "family_id" in sanity.columns else set()
    if norm.empty or "family_cds_qc_passed" not in norm.columns or not keep_families:
        pd.DataFrame(columns=manifest_columns).to_csv(manifest_out, sep="\t", index=False)
        log("Stage 9 FASTA writing complete: fastas=0")
        return
    norm = norm[
        (norm["family_cds_qc_passed"] == True)
        & (norm["passed"] == True)
        & (norm["family_id"].isin(keep_families))
    ].copy()
    order = [s.token for s in species]
    order_rank = {token: i for i, token in enumerate(order)}
    expected_species_count = len(species)
    manifest_rows: List[Dict[str, Any]] = []
    used_names: Dict[str, str] = {}

    for family_id, fam in norm.groupby("family_id", sort=False):
        fam = fam.copy()
        mode = str(fam["orthology_mode"].iloc[0]) if "orthology_mode" in fam.columns else "strict_singleton"
        if mode != "reference_gene_1to1_present_species" and fam["token"].nunique() != expected_species_count:
            continue
        reference_symbol = str(
            fam["reference_symbol"].iloc[0]
            if "reference_symbol" in fam.columns
            else fam["human_symbol"].iloc[0]
        )
        reference_geneid = int(
            fam["reference_GeneID"].iloc[0]
            if "reference_GeneID" in fam.columns
            else fam["human_GeneID"].iloc[0]
        )
        reference_taxid = int(fam["reference_taxid"].iloc[0]) if "reference_taxid" in fam.columns else 9606
        reference_token = str(fam["reference_token"].iloc[0]) if "reference_token" in fam.columns else "human"
        base = safe_symbol(reference_symbol)
        filename = f"{base}.reference_1to1.cds.fasta" if mode == "reference_gene_1to1_present_species" else f"{base}.fasta"
        if filename in used_names and used_names[filename] != family_id:
            if mode == "reference_gene_1to1_present_species":
                filename = f"{base}__{reference_geneid}.reference_1to1.cds.fasta"
            else:
                filename = f"{base}__{reference_geneid}.fasta"
        used_names[filename] = family_id
        if mode == "reference_gene_1to1_present_species":
            meta_filename = filename.replace(".reference_1to1.cds.fasta", ".reference_1to1.meta.tsv")
            rejected_filename = filename.replace(".reference_1to1.cds.fasta", ".reference_1to1.rejected.tsv")
        else:
            meta_filename = filename.replace(".fasta", ".meta.tsv")
            rejected_filename = filename.replace(".fasta", ".rejected.tsv")
        fasta_path = FASTAS / filename
        meta_path = FASTAS / meta_filename
        fam["rank"] = fam["token"].map(order_rank)
        fam = fam.sort_values("rank")
        with fasta_path.open("w") as fh:
            for row in fam.to_dict("records"):
                fh.write(f">{row['token']}\n")
                fh.write(wrap_sequence(str(row["normalized_sequence"])) + "\n")

        meta_rows: List[Dict[str, Any]] = []
        for row in fam.to_dict("records"):
            meta = assembly_meta.get(row["token"], {})
            meta_rows.append(
                {
                    "family_id": family_id,
                    "orthology_mode": mode,
                    "reference_symbol": reference_symbol,
                    "reference_GeneID": reference_geneid,
                    "reference_taxid": reference_taxid,
                    "reference_token": reference_token,
                    "human_symbol": reference_symbol,
                    "fasta_path": rel_path(fasta_path),
                    "token": row["token"],
                    "taxid": row["taxid"],
                    "GeneID": row["GeneID"],
                    "species_symbol": row.get("species_symbol", ""),
                    "transcript_accession": row.get("transcript_accession", ""),
                    "cds_accession": row.get("cds_accession", ""),
                    "protein_accession": row.get("protein_accession", ""),
                    "cds_length_original": row.get("cds_length_original", ""),
                    "cds_length_normalized": row.get("cds_length_normalized", ""),
                    "protein_length": row.get("protein_length", ""),
                    "select_category": row.get("select_category", ""),
                    "accession_class": row.get("accession_class", ""),
                    "had_terminal_stop": row.get("had_terminal_stop", ""),
                    "terminal_stop_removed": row.get("terminal_stop_removed", ""),
                    "gcf_accession": meta.get("gcf_accession", ""),
                    "annotation_release_id": meta.get("annotation_release_id", ""),
                }
            )
        pd.DataFrame(meta_rows).to_csv(meta_path, sep="\t", index=False)
        if mode == "reference_gene_1to1_present_species":
            query_dir = reference_query_dir(reference_symbol, reference_geneid)
            query_dir.mkdir(parents=True, exist_ok=True)
            result_fasta = query_dir / filename
            result_meta = query_dir / meta_filename
            shutil.copyfile(fasta_path, result_fasta)
            shutil.copyfile(meta_path, result_meta)
            rejected_path = ORTHOLOGY / "rejected.tsv"
            if rejected_path.exists():
                shutil.copyfile(rejected_path, query_dir / rejected_filename)
        manifest_rows.append(
            {
                "family_id": family_id,
                "orthology_mode": mode,
                "reference_symbol": reference_symbol,
                "reference_GeneID": reference_geneid,
                "reference_taxid": reference_taxid,
                "reference_token": reference_token,
                "human_symbol": reference_symbol,
                "human_GeneID": reference_geneid,
                "fasta_path": rel_path(fasta_path),
                "meta_path": rel_path(meta_path),
                "sequence_count": len(fam),
                "sanity_status": sanity.loc[sanity["family_id"] == family_id, "status"].iloc[0],
            }
        )

    manifest_df = pd.DataFrame(manifest_rows, columns=manifest_columns)
    manifest_df.to_csv(manifest_out, sep="\t", index=False)
    log(f"Stage 9 FASTA writing complete: fastas={len(manifest_df):,}")


def stage10_report() -> None:
    summary: Dict[str, Any] = {"generated_at": datetime.now().isoformat()}
    paths = {
        "assembly_index": INDEXES / "assembly_index.parquet",
        "gene_index": INDEXES / "gene_index.parquet",
        "cds_index": INDEXES / "cds_index.parquet",
        "ortholog_edges": ORTHOLOGY / "ortholog_edges.parquet",
        "strict_singleton": ORTHOLOGY / "strict_singleton.parquet",
        "representative_cds": SELECTION / "representative_cds.parquet",
        "cds_qc": QC / "cds_qc.parquet",
        "family_sanity": QC / "family_sanity.parquet",
    }
    for name, path in paths.items():
        if path.exists():
            df = pd.read_parquet(path)
            summary[name] = {"rows": int(len(df))}
            for col in ["family_id", "taxid", "GeneID"]:
                if col in df.columns:
                    summary[name][f"unique_{col}"] = int(df[col].nunique())
    manifest = FASTAS / "manifest.tsv"
    mf = pd.DataFrame()
    if manifest.exists():
        mf = pd.read_csv(manifest, sep="\t")
        summary["fastas"] = {"count": int(len(mf))}
    rejected = ORTHOLOGY / "rejected.parquet"
    if rejected.exists():
        rej = pd.read_parquet(rejected)
        reason_col = "reason" if "reason" in rej.columns else "reject_reason"
        summary["rejected_reasons"] = rej[reason_col].value_counts().to_dict() if (not rej.empty and reason_col in rej.columns) else {}
    qc_path = QC / "cds_qc.parquet"
    if qc_path.exists():
        qc_df = pd.read_parquet(qc_path)
        summary["cds_qc_fail_reasons"] = (
            qc_df.loc[~qc_df["passed"], "fail_reason"].value_counts().to_dict()
            if not qc_df.empty
            else {}
        )
    write_json(REPORTS / "summary.json", summary)

    if not mf.empty and "orthology_mode" in mf.columns:
        ref_mode = mf[mf["orthology_mode"] == "reference_gene_1to1_present_species"]
        if not ref_mode.empty:
            selected = pd.read_parquet(SELECTION / "representative_cds.parquet") if (SELECTION / "representative_cds.parquet").exists() else pd.DataFrame()
            norm = pd.read_parquet(SELECTION / "representative_cds.normalized.parquet") if (SELECTION / "representative_cds.normalized.parquet").exists() else pd.DataFrame()
            qc = pd.read_parquet(QC / "cds_qc.parquet") if (QC / "cds_qc.parquet").exists() else pd.DataFrame()
            rejected = pd.read_parquet(ORTHOLOGY / "rejected.parquet") if (ORTHOLOGY / "rejected.parquet").exists() else pd.DataFrame()
            for row in ref_mode.to_dict("records"):
                query_dir = reference_query_dir(str(row["reference_symbol"]), int(row["reference_GeneID"]))
                query_dir.mkdir(parents=True, exist_ok=True)
                family_id = row["family_id"]
                if not selected.empty:
                    selected[selected["family_id"] == family_id].to_csv(query_dir / "selected_cds.tsv", sep="\t", index=False)
                if not norm.empty:
                    norm_family = norm[norm["family_id"] == family_id]
                    norm_family.to_csv(query_dir / "cds_qc_input.tsv", sep="\t", index=False)
                    norm_family[norm_family["passed"] == True].to_csv(query_dir / "cds_pass.tsv", sep="\t", index=False)
                    norm_family[norm_family["passed"] != True].to_csv(query_dir / "cds_rejected.tsv", sep="\t", index=False)
                if not qc.empty:
                    qc[qc["family_id"] == family_id].to_csv(query_dir / "cds_qc.tsv", sep="\t", index=False)
                if not rejected.empty:
                    rejected.to_csv(query_dir / f"{row['reference_symbol']}.reference_1to1.rejected.tsv", sep="\t", index=False)
                query_summary = {
                    "generated_at": summary["generated_at"],
                    "orthology_mode": "reference_gene_1to1_present_species",
                    "reference_symbol": row["reference_symbol"],
                    "reference_GeneID": int(row["reference_GeneID"]),
                    "reference_taxid": int(row["reference_taxid"]),
                    "manifest_species_count": int(len(read_manifest())),
                    "final_fasta_sequence_count": int(row["sequence_count"]),
                    "fasta_path": row["fasta_path"],
                    "meta_path": row["meta_path"],
                    "rejected_reasons": rejected["reject_reason"].value_counts().to_dict()
                    if (not rejected.empty and "reject_reason" in rejected.columns)
                    else {},
                    "cds_qc_fail_reasons": qc.loc[~qc["passed"], "fail_reason"].value_counts().to_dict()
                    if (not qc.empty and "passed" in qc.columns)
                    else {},
                }
                write_json(query_dir / "summary.json", query_summary)
                (query_dir / "summary.html").write_text(
                    "<html><head><meta charset='utf-8'><title>Reference Gene 1:1 Summary</title></head><body>"
                    f"<h1>{row['reference_symbol']} reference-gene 1:1 summary</h1>"
                    f"<p>Final FASTA sequences: {row['sequence_count']}</p>"
                    f"<p>FASTA: {row['fasta_path']}</p>"
                    "</body></html>"
                )

    html = [
        "<html><head><meta charset='utf-8'><title>NCBI Ortholog CDS Summary</title>",
        "<style>body{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;max-width:1100px;margin:40px auto;line-height:1.45}table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:4px 8px;text-align:left}</style>",
        "</head><body>",
        "<h1>NCBI Ortholog CDS Summary</h1>",
        f"<p>Generated at {summary['generated_at']}</p>",
        "<h2>Tables</h2><table><tr><th>Name</th><th>Rows</th><th>Unique families</th></tr>",
    ]
    for name, data in summary.items():
        if isinstance(data, dict) and "rows" in data:
            html.append(
                f"<tr><td>{name}</td><td>{data.get('rows','')}</td><td>{data.get('unique_family_id','')}</td></tr>"
            )
    html.append("</table>")
    if "fastas" in summary:
        html.append(f"<h2>FASTAs</h2><p>{summary['fastas']['count']} files generated.</p>")
    if "rejected_reasons" in summary:
        html.append("<h2>Rejected Components</h2><table><tr><th>Reason</th><th>Count</th></tr>")
        for reason, count in summary["rejected_reasons"].items():
            html.append(f"<tr><td>{reason}</td><td>{count}</td></tr>")
        html.append("</table>")
    html.append("</body></html>")
    (REPORTS / "summary.html").write_text("\n".join(html))
    log("Stage 10 summary report complete")


def parse_steps(raw: str) -> List[str]:
    if raw == "all":
        return [
            "preflight",
            "bulk",
            "assemblies",
            "indexes",
            "edges",
            "singletons",
            "select",
            "qc",
            "sanity",
            "fastas",
            "report",
        ]
    return [x.strip() for x in raw.split(",") if x.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--steps",
        default="all",
        help="Comma-separated steps or all. Steps: preflight,bulk,assemblies,indexes,edges,singletons,select,qc,sanity,fastas,report",
    )
    parser.add_argument(
        "--reference-taxid",
        type=int,
        default=9606,
        help="TaxID whose gene symbol is used for family IDs and FASTA filenames (default: 9606, human).",
    )
    parser.add_argument(
        "--orthology-mode",
        choices=["strict_singleton", "reference_gene_1to1_present_species"],
        default="strict_singleton",
        help="Orthology parsing mode. strict_singleton keeps only all-species N-way singletons; "
        "reference_gene_1to1_present_species keeps species with clear 1:1 orthologs to one query gene.",
    )
    parser.add_argument("--reference-symbol", help="Reference-species gene symbol for reference_gene_1to1_present_species mode")
    parser.add_argument("--reference-gene-id", type=int, help="Reference-species NCBI GeneID for reference_gene_1to1_present_species mode")
    parser.add_argument("--min-sequences", type=int, default=2, help="Minimum sequences required for reference_gene_1to1_present_species output")
    parser.add_argument(
        "--exclude-reference",
        action="store_true",
        help="Do not include the reference species sequence in reference_gene_1to1_present_species FASTA output",
    )
    parser.add_argument("--manifest", default=str(MANIFEST), help="Species manifest TSV")
    parser.add_argument("--input-root", help="Input root containing ncbi_bulk/ and assembly_packages/")
    parser.add_argument("--output-root", help="Output root for indexes, reports, FASTA files and QC")
    parser.add_argument("--datasets", default=str(DATASETS), help="NCBI Datasets CLI executable")
    parser.add_argument("--offline", action="store_true", help="Use local input-root fixtures; do not download NCBI data")
    parser.add_argument("--force", action="store_true", help="Rebuild selected stages")
    args = parser.parse_args()

    configure_paths(
        manifest=args.manifest,
        input_root=args.input_root,
        output_root=args.output_root,
        datasets=args.datasets,
    )
    os.chdir(PROJECT_ROOT)
    ensure_dirs()
    species = read_manifest()
    manifest_taxids = {s.taxid for s in species}
    if args.reference_taxid not in manifest_taxids:
        raise ValueError(
            f"--reference-taxid {args.reference_taxid} is not present in config/species_manifest.tsv"
        )
    steps = parse_steps(args.steps)
    orth_path = latest_bulk_file("gene_orthologs")
    info_path = latest_bulk_file("gene_info")

    if "preflight" in steps:
        stage0_preflight(species, args.force, args.offline)
    if "bulk" in steps:
        orth_path, info_path = stage1_bulk(args.force, args.offline)
    needs_bulk_downstream = any(
        step in steps
        for step in ["indexes", "edges", "singletons", "select", "qc", "sanity", "fastas", "report"]
    )
    if needs_bulk_downstream and (not orth_path or not info_path):
        orth_path, info_path = stage1_bulk(False, args.offline)
    assembly_species = species
    if (
        args.orthology_mode == "reference_gene_1to1_present_species"
        and any(step in steps for step in ["assemblies", "indexes", "singletons", "select", "qc", "sanity", "fastas", "report"])
    ):
        if not orth_path or not info_path:
            raise RuntimeError("Reference-gene mode requires NCBI Gene orthology and gene_info bulk files")
        assembly_species = species_needed_for_reference_gene_mode(
            species,
            orth_path,
            info_path,
            args.reference_taxid,
            args.reference_symbol,
            args.reference_gene_id,
        )
    assembly_include = (
        "cds,gff3,seq-report"
        if args.orthology_mode == "reference_gene_1to1_present_species"
        else "cds,protein,rna,gff3,seq-report"
    )
    if "assemblies" in steps:
        stage2_download_assemblies(assembly_species, args.force, args.offline, assembly_include)
    if "indexes" in steps:
        stage3_build_indexes(assembly_species, info_path, args.force)
    if "edges" in steps:
        stage4_orthology_edges(
            species,
            orth_path,
            args.force,
            args.reference_taxid
            if args.orthology_mode == "reference_gene_1to1_present_species"
            else None,
        )
    if "singletons" in steps:
        if args.orthology_mode == "reference_gene_1to1_present_species":
            stage5_reference_gene_1to1_present_species(
                species,
                args.reference_taxid,
                args.reference_symbol,
                args.reference_gene_id,
                not args.exclude_reference,
                args.min_sequences,
                args.force,
            )
        else:
            stage5_strict_singletons(species, args.reference_taxid, args.force)
    if "select" in steps:
        stage6_select_representatives(len(species), args.force)
    if "qc" in steps:
        stage7_cds_qc(len(species), args.force)
    if "sanity" in steps:
        stage8_family_sanity(args.force)
    if "fastas" in steps:
        stage9_write_fastas(species, args.force)
    if "report" in steps:
        stage10_report()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        raise
    except Exception as exc:
        log(f"ERROR: {exc}")
        raise
