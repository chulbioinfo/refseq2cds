#!/usr/bin/env python3
"""
Tree-free codon alignments with MAFFT + PAL2NAL.

This driver takes refseq2cds FASTA outputs, translates each CDS to amino acids,
aligns proteins with MAFFT, and converts the protein MSA plus original CDS
sequences into a codon alignment with PAL2NAL.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from Bio import SeqIO
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.Seq import Seq


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_DIR = ROOT / "fastas"
DEFAULT_OUTPUT_DIR = ROOT / "alignments" / "mafft_pal2nal"
DEFAULT_FASTA_MANIFEST = ROOT / "fastas" / "manifest.tsv"
DEFAULT_SPECIES_MANIFEST = ROOT / "config" / "species_manifest.tsv"

STANDARD_TABLE = unambiguous_dna_by_id[1]
STOP_CODONS = set(STANDARD_TABLE.stop_codons)


@dataclass(frozen=True)
class AlignmentResult:
    status: str
    family_id: str
    symbol: str
    input_fasta: str
    codon_alignment_path: str = ""
    aa_fasta_path: str = ""
    aa_alignment_path: str = ""
    log_path: str = ""
    map_paths: str = ""
    sequence_count: int = 0
    alignment_length_nt: int = 0
    alignment_length_codons: int = 0
    max_gap_fraction: float = 0.0
    reason: str = ""


def log(message: str) -> None:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {message}", flush=True)


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def wrap_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
    for record in SeqIO.parse(str(path), "fasta"):
        ident = record.id.split()[0]
        if ident in records:
            raise ValueError(f"Duplicate FASTA ID {ident} in {path}")
        records[ident] = str(record.seq).upper()
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def write_fasta(records: Dict[str, str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for ident, seq in records.items():
            fh.write(f">{ident}\n")
            fh.write(wrap_sequence(seq) + "\n")


def read_tsv(path: Path) -> List[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def expected_tokens_from_manifest(path: Path) -> List[str]:
    rows = read_tsv(path)
    return [row["token"] for row in rows if row.get("token")]


def family_rows_from_manifest(path: Path) -> Dict[str, dict]:
    rows = read_tsv(path)
    out: Dict[str, dict] = {}
    for row in rows:
        fasta_path = row.get("fasta_path", "")
        if fasta_path:
            stem = Path(fasta_path).stem
            out[stem] = row
            if row.get("reference_symbol"):
                out[row["reference_symbol"]] = row
            if row.get("human_symbol"):
                out[row["human_symbol"]] = row
    return out


def symbol_from_fasta_stem(stem: str) -> str:
    for suffix in [".reference_1to1.cds", ".cds"]:
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def flatten_symbols(values: Optional[Sequence[str]]) -> List[str]:
    if not values:
        return []
    symbols: List[str] = []
    for value in values:
        for part in value.split(","):
            part = part.strip()
            if part:
                symbols.append(symbol_from_fasta_stem(Path(part).stem))
    return symbols


def collect_fastas(input_dir: Path, symbols: Optional[Sequence[str]], limit: Optional[int]) -> List[Path]:
    candidates = sorted(
        p
        for p in input_dir.glob("*.fasta")
        if p.is_file() and not p.name.startswith("._")
    )
    wanted = set(flatten_symbols(symbols))
    if wanted:
        candidates = [p for p in candidates if p.stem in wanted or symbol_from_fasta_stem(p.stem) in wanted]
        observed = {p.stem for p in candidates} | {symbol_from_fasta_stem(p.stem) for p in candidates}
        missing = sorted(wanted - observed)
        if missing:
            raise FileNotFoundError(f"Requested symbols not found in {input_dir}: {missing[:10]}")
    if limit is not None:
        candidates = candidates[:limit]
    if not candidates:
        raise FileNotFoundError(f"No FASTA files selected from {input_dir}")
    return candidates


def resolve_executable(name: str, label: str) -> str:
    path = Path(name)
    if path.exists():
        return str(path)
    found = shutil.which(name)
    if found:
        return found
    raise FileNotFoundError(
        f"Missing {label}: {name}. Install it or pass an explicit path with the corresponding option."
    )


def pal2nal_command(pal2nal: str) -> List[str]:
    resolved = resolve_executable(pal2nal, "PAL2NAL")
    path = Path(resolved)
    if path.suffix == ".pl" and not os.access(path, os.X_OK):
        perl = shutil.which("perl")
        if not perl:
            raise FileNotFoundError("PAL2NAL is a Perl script, but perl was not found on PATH")
        return [perl, resolved]
    return [resolved]


def command_text(cmd: Sequence[str]) -> str:
    return " ".join(str(x) for x in cmd)


def command_version(cmd: Sequence[str]) -> str:
    probes = [list(cmd) + ["--version"], list(cmd) + ["-h"]]
    for probe in probes:
        try:
            proc = subprocess.run(probe, capture_output=True, text=True, timeout=10)
        except Exception:
            continue
        text = (proc.stdout + "\n" + proc.stderr).strip()
        if text:
            return text.splitlines()[0][:200]
    return "unknown"


def translate_cds_records(nuc_records: Dict[str, str]) -> Dict[str, str]:
    aa_records: Dict[str, str] = {}
    for ident, seq in nuc_records.items():
        clean = seq.upper().replace("U", "T").replace("-", "")
        if len(clean) % 3 != 0:
            raise ValueError(f"{ident} CDS length is not a multiple of 3: {len(clean)}")
        aa = str(Seq(clean).translate(table=1, to_stop=False))
        if "*" in aa:
            raise ValueError(f"{ident} translates with stop codon(s); run CDS QC before alignment")
        aa_records[ident] = aa
    return aa_records


def run_mafft(mafft: str, threads: int, aa_fasta: Path, aa_alignment: Path, log_path: Path) -> None:
    cmd = [mafft, "--auto", "--thread", str(threads), str(aa_fasta)]
    with aa_alignment.open("w") as out, log_path.open("a") as err:
        err.write(f"\n# MAFFT\n{command_text(cmd)}\n")
        proc = subprocess.run(cmd, stdout=out, stderr=err, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"MAFFT failed with exit code {proc.returncode}")


def run_pal2nal(
    pal2nal_cmd: Sequence[str],
    aa_alignment: Path,
    nuc_fasta: Path,
    codon_alignment: Path,
    log_path: Path,
    codon_table: str,
) -> None:
    codon_table_arg = {"universal": "1", "vmitochondria": "2"}[codon_table]
    cmd = list(pal2nal_cmd) + [
        str(aa_alignment),
        str(nuc_fasta),
        "-output",
        "fasta",
        "-codontable",
        codon_table_arg,
    ]
    tmp_out = codon_alignment.with_suffix(codon_alignment.suffix + ".tmp")
    with tmp_out.open("w") as out, log_path.open("a") as err:
        err.write(f"\n# PAL2NAL\n{command_text(cmd)}\n")
        proc = subprocess.run(cmd, stdout=out, stderr=err, text=True)
    if proc.returncode != 0:
        tmp_out.unlink(missing_ok=True)
        raise RuntimeError(f"PAL2NAL failed with exit code {proc.returncode}")
    tmp_out.replace(codon_alignment)


def validate_codon_alignment(codon_alignment: Path, expected_ids: Iterable[str]) -> Dict[str, object]:
    records = read_fasta(codon_alignment)
    expected = set(expected_ids)
    observed = set(records)
    if observed != expected:
        raise ValueError(f"Alignment IDs do not match input IDs: {sorted(observed ^ expected)}")
    lengths = {len(seq) for seq in records.values()}
    if len(lengths) != 1:
        raise ValueError(f"Alignment sequences have different lengths: {sorted(lengths)}")
    alignment_length = lengths.pop()
    if alignment_length % 3 != 0:
        raise ValueError(f"Codon alignment length is not a multiple of 3: {alignment_length}")
    partial_gap_codons = 0
    stop_codons = 0
    max_gap_fraction = 0.0
    for ident, seq in records.items():
        seq = seq.upper()
        gap_fraction = seq.count("-") / len(seq) if seq else 0.0
        max_gap_fraction = max(max_gap_fraction, gap_fraction)
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            if codon == "---":
                continue
            if "-" in codon:
                partial_gap_codons += 1
                continue
            if "N" in codon:
                continue
            if codon in STOP_CODONS:
                stop_codons += 1
    if partial_gap_codons:
        raise ValueError(f"Alignment contains partial codon gaps: {partial_gap_codons}")
    if stop_codons:
        raise ValueError(f"Alignment contains stop codons: {stop_codons}")
    return {
        "sequence_count": len(records),
        "alignment_length_nt": alignment_length,
        "alignment_length_codons": alignment_length // 3,
        "max_gap_fraction": max_gap_fraction,
    }


def write_token_map(codon_alignment: Path, family_id: str, token: str, out_path: Path) -> Optional[Path]:
    records = read_fasta(codon_alignment)
    if token not in records:
        return None
    seq = records[token].upper()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cds_codon_index = 0
    with gzip.open(out_path, "wt", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "family_id",
                "token",
                "aln_codon_index_1based",
                "cds_codon_index_1based",
                "cds_nt_start_1based",
                "cds_nt_end_1based",
                "is_gap",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for aln_codon_index, i in enumerate(range(0, len(seq), 3), start=1):
            codon = seq[i : i + 3]
            is_gap = codon == "---"
            if is_gap:
                cds_index = ""
                nt_start = ""
                nt_end = ""
            else:
                cds_codon_index += 1
                cds_index = cds_codon_index
                nt_start = (cds_codon_index - 1) * 3 + 1
                nt_end = cds_codon_index * 3
            writer.writerow(
                {
                    "family_id": family_id,
                    "token": token,
                    "aln_codon_index_1based": aln_codon_index,
                    "cds_codon_index_1based": cds_index,
                    "cds_nt_start_1based": nt_start,
                    "cds_nt_end_1based": nt_end,
                    "is_gap": str(is_gap).lower(),
                }
            )
    return out_path


def process_one(
    fasta_path: Path,
    output_dir: Path,
    mafft: str,
    pal2nal_cmd: Sequence[str],
    threads_per_mafft: int,
    codon_table: str,
    expected_tokens: Sequence[str],
    family_manifest: Dict[str, dict],
    map_tokens: Sequence[str],
    force: bool,
) -> AlignmentResult:
    stem = fasta_path.stem
    manifest_row = family_manifest.get(stem, {})
    family_id = manifest_row.get("family_id", stem)
    symbol = manifest_row.get("reference_symbol") or manifest_row.get("human_symbol") or stem
    codon_dir = output_dir / "codon"
    aa_dir = output_dir / "aa"
    aa_aln_dir = output_dir / "aa_aligned"
    log_dir = output_dir / "logs"
    map_dir = output_dir / "maps"
    codon_alignment = codon_dir / f"{stem}.codon.fasta"
    aa_fasta = aa_dir / f"{stem}.aa.fasta"
    aa_alignment = aa_aln_dir / f"{stem}.aa.aln.fasta"
    log_path = log_dir / f"{stem}.log"
    for directory in [codon_dir, aa_dir, aa_aln_dir, log_dir, map_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    try:
        if codon_alignment.exists() and not force:
            qc = validate_codon_alignment(codon_alignment, read_fasta(fasta_path).keys())
        else:
            log_path.write_text(f"# refseq2cds mafft-pal2nal alignment for {stem}\n")
            nuc_records = read_fasta(fasta_path)
            if expected_tokens:
                observed_tokens = set(nuc_records)
                expected_token_set = set(expected_tokens)
                unknown_tokens = observed_tokens - expected_token_set
                if unknown_tokens:
                    raise ValueError(f"{fasta_path} IDs are not in the species manifest: {sorted(unknown_tokens)}")
                mode = manifest_row.get("orthology_mode", "strict_singleton")
                if mode != "reference_gene_1to1_present_species" and observed_tokens != expected_token_set:
                    raise ValueError(
                        f"{fasta_path} IDs do not match manifest tokens: "
                        f"{sorted(observed_tokens ^ expected_token_set)}"
                    )
            aa_records = translate_cds_records(nuc_records)
            write_fasta(aa_records, aa_fasta)
            run_mafft(mafft, threads_per_mafft, aa_fasta, aa_alignment, log_path)
            run_pal2nal(pal2nal_cmd, aa_alignment, fasta_path, codon_alignment, log_path, codon_table)
            qc = validate_codon_alignment(codon_alignment, nuc_records.keys())

        map_paths: List[str] = []
        if map_tokens:
            for token in map_tokens:
                map_path = map_dir / f"{stem}.{token}.codon_map.tsv.gz"
                if map_path.exists() and not force:
                    map_paths.append(rel(map_path))
                    continue
                written = write_token_map(codon_alignment, family_id, token, map_path)
                if written is not None:
                    map_paths.append(rel(written))

        return AlignmentResult(
            status="pass",
            family_id=family_id,
            symbol=symbol,
            input_fasta=rel(fasta_path),
            codon_alignment_path=rel(codon_alignment),
            aa_fasta_path=rel(aa_fasta),
            aa_alignment_path=rel(aa_alignment),
            log_path=rel(log_path),
            map_paths=";".join(map_paths),
            sequence_count=int(qc["sequence_count"]),
            alignment_length_nt=int(qc["alignment_length_nt"]),
            alignment_length_codons=int(qc["alignment_length_codons"]),
            max_gap_fraction=float(qc["max_gap_fraction"]),
        )
    except Exception as exc:
        with log_path.open("a") as fh:
            fh.write(f"\n# FAILED\n{type(exc).__name__}: {exc}\n")
        return AlignmentResult(
            status="fail",
            family_id=family_id,
            symbol=symbol,
            input_fasta=rel(fasta_path),
            codon_alignment_path=rel(codon_alignment),
            aa_fasta_path=rel(aa_fasta),
            aa_alignment_path=rel(aa_alignment),
            log_path=rel(log_path),
            reason=f"{type(exc).__name__}: {exc}",
        )


def write_table(path: Path, rows: List[dict], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create tree-free codon alignments from refseq2cds FASTA files using MAFFT + PAL2NAL.",
    )
    parser.add_argument("--input-dir", default=str(DEFAULT_INPUT_DIR), help="Directory containing refseq2cds FASTA files")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="Alignment output directory")
    parser.add_argument("--fasta-manifest", default=str(DEFAULT_FASTA_MANIFEST), help="fastas/manifest.tsv path")
    parser.add_argument("--species-manifest", default=str(DEFAULT_SPECIES_MANIFEST), help="config/species_manifest.tsv path")
    parser.add_argument("--mafft", default="mafft", help="MAFFT executable path or command name")
    parser.add_argument("--pal2nal", default="pal2nal.pl", help="PAL2NAL executable path or command name")
    parser.add_argument("--threads-per-mafft", type=int, default=1, help="Threads passed to each MAFFT process")
    parser.add_argument("--jobs", type=int, default=1, help="Number of families to align concurrently")
    parser.add_argument("--limit", type=int, help="Only process the first N selected FASTA files")
    parser.add_argument("--symbols", nargs="*", help="Only process selected symbols/stems; comma-separated values are accepted")
    parser.add_argument("--map-token", action="append", default=[], help="Write alignment-to-CDS codon map for this token")
    parser.add_argument("--codon-table", choices=["universal", "vmitochondria"], default="universal")
    parser.add_argument("--force", action="store_true", help="Overwrite existing alignments")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_parser().parse_args(argv)
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    fasta_manifest_path = Path(args.fasta_manifest)
    species_manifest_path = Path(args.species_manifest)
    mafft = resolve_executable(args.mafft, "MAFFT")
    pal2nal_cmd = pal2nal_command(args.pal2nal)
    expected_tokens = expected_tokens_from_manifest(species_manifest_path)
    family_manifest = family_rows_from_manifest(fasta_manifest_path)
    fastas = collect_fastas(input_dir, args.symbols, args.limit)
    output_dir.mkdir(parents=True, exist_ok=True)
    log(
        "Starting MAFFT+PAL2NAL alignments: "
        f"fastas={len(fastas):,}; jobs={args.jobs}; threads_per_mafft={args.threads_per_mafft}"
    )

    results: List[AlignmentResult] = []
    if args.jobs == 1:
        for fasta_path in fastas:
            results.append(
                process_one(
                    fasta_path,
                    output_dir,
                    mafft,
                    pal2nal_cmd,
                    args.threads_per_mafft,
                    args.codon_table,
                    expected_tokens,
                    family_manifest,
                    args.map_token,
                    args.force,
                )
            )
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [
                pool.submit(
                    process_one,
                    fasta_path,
                    output_dir,
                    mafft,
                    pal2nal_cmd,
                    args.threads_per_mafft,
                    args.codon_table,
                    expected_tokens,
                    family_manifest,
                    args.map_token,
                    args.force,
                )
                for fasta_path in fastas
            ]
            for future in as_completed(futures):
                results.append(future.result())

    rows = [result.__dict__ for result in sorted(results, key=lambda r: r.input_fasta)]
    pass_rows = [row for row in rows if row["status"] == "pass"]
    fail_rows = [row for row in rows if row["status"] != "pass"]
    fieldnames = list(AlignmentResult.__dataclass_fields__.keys())
    write_table(output_dir / "manifest.tsv", pass_rows, fieldnames)
    write_table(output_dir / "failed.tsv", fail_rows, fieldnames)
    summary = {
        "generated_at": datetime.now().isoformat(),
        "mode": "mafft-pal2nal",
        "input_dir": rel(input_dir),
        "output_dir": rel(output_dir),
        "selected_fastas": len(fastas),
        "passed": len(pass_rows),
        "failed": len(fail_rows),
        "mafft": mafft,
        "pal2nal": command_text(pal2nal_cmd),
        "mafft_version": command_version([mafft]),
        "pal2nal_version": command_version(list(pal2nal_cmd)),
        "threads_per_mafft": args.threads_per_mafft,
        "jobs": args.jobs,
        "codon_table": args.codon_table,
        "map_tokens": args.map_token,
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    log(f"Done: passed={len(pass_rows):,}; failed={len(fail_rows):,}; output={output_dir}")
    if fail_rows:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
