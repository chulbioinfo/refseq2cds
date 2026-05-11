#!/usr/bin/env python3
"""
refseq2cds

Command-line facade for the assembly-exact NCBI singleton CDS pipeline.

This module keeps the project usable as both:
1. a repository-local workflow runner, and
2. an installable Python console script.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional, Sequence


MODULE_ROOT = Path(__file__).resolve().parent


def find_project_root() -> Path:
    """Find the repository root containing workflow/scripts.

    Editable installs resolve to the repository path. Non-editable installs can
    still be used from a cloned checkout because the current working directory
    is checked first.
    """
    candidates = [Path.cwd().resolve(), MODULE_ROOT]
    for candidate in list(candidates):
        candidates.extend(candidate.parents)
    seen = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        if (candidate / "workflow" / "scripts" / "run_cds_pipeline.py").exists():
            return candidate
    return MODULE_ROOT


ROOT = find_project_root()
DEFAULT_MANIFEST = ROOT / "config" / "species_manifest.tsv"
RUN_DRIVER = ROOT / "workflow" / "scripts" / "run_cds_pipeline.py"
MATRIX_DRIVER = ROOT / "workflow" / "scripts" / "build_human_cds_position_matrices.py"
ALIGN_DRIVER = ROOT / "workflow" / "scripts" / "align_mafft_pal2nal.py"

KNOWN_TAXID_LABELS = {
    9606: ("human", "human"),
    9598: ("chimpanzee", "chimpanzee"),
    9597: ("bonobo", "bonobo"),
    9595: ("gorilla", "gorilla"),
    9601: ("Sumatran_orangutan", "Sumatran orangutan"),
    9600: ("Bornean_orangutan", "Bornean orangutan"),
    9590: ("siamang_gibbon", "siamang gibbon"),
    9541: ("crab-eating_macaque", "crab-eating macaque"),
    9545: ("pig-tailed_macaque", "pig-tailed macaque"),
    9483: ("common_marmoset", "common marmoset"),
    27679: ("Bolivian_squirrel_monkey", "Bolivian squirrel monkey"),
    9470: ("sunda_slow_loris", "sunda slow loris"),
    9447: ("ring-tailed_lemur", "ring-tailed lemur"),
    110931: ("Philippine_flying_lemur", "Philippine flying lemur"),
}


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def run(cmd: Sequence[str]) -> None:
    log("RUN " + " ".join(cmd))
    subprocess.run(cmd, cwd=ROOT, check=True)


def run_capture(cmd: Sequence[str]) -> str:
    log("RUN " + " ".join(cmd))
    proc = subprocess.run(cmd, cwd=ROOT, check=True, capture_output=True, text=True)
    return proc.stdout


def python_executable() -> str:
    return sys.executable


def require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")


def download_datasets_cli() -> None:
    bin_dir = ROOT / "bin"
    bin_dir.mkdir(exist_ok=True)
    system = platform.system().lower()
    if system == "darwin":
        os_name = "mac"
    elif system == "linux":
        os_name = "linux"
    else:
        raise RuntimeError(f"Unsupported OS for automatic NCBI Datasets download: {platform.system()}")

    for tool in ["datasets", "dataformat"]:
        out = bin_dir / tool
        if out.exists() and os.access(out, os.X_OK):
            log(f"{tool} already exists: {out}")
            continue
        url = f"https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/{os_name}/{tool}"
        run(["curl", "-L", "--fail", "-o", str(out), url])
        out.chmod(0o755)
        log(f"Installed {tool}: {out}")


def cmd_init_manifest(args: argparse.Namespace) -> None:
    target = Path(args.output)
    if target.exists() and not args.force:
        raise FileExistsError(f"{target} already exists. Use --force to overwrite.")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        "\n".join(
            [
                "token\ttaxid\tscientific_name\tcommon_name\tgcf_accession\toutgroup",
                "human\t9606\tHomo sapiens\thuman\tGCF_009914755.1\tfalse",
                "chimpanzee\t9598\tPan troglodytes\tchimpanzee\tGCF_028858775.2\tfalse",
                "bonobo\t9597\tPan paniscus\tbonobo\tGCF_029289425.2\tfalse",
                "gorilla\t9595\tGorilla gorilla\tgorilla\tGCF_029281585.2\tfalse",
                "Sumatran_orangutan\t9601\tPongo abelii\tSumatran orangutan\tGCF_028885655.2\tfalse",
                "Bornean_orangutan\t9600\tPongo pygmaeus\tBornean orangutan\tGCF_028885625.2\tfalse",
                "siamang_gibbon\t9590\tSymphalangus syndactylus\tsiamang gibbon\tGCF_028878055.3\tfalse",
                "crab-eating_macaque\t9541\tMacaca fascicularis\tcrab-eating macaque\tGCF_037993035.2\tfalse",
                "pig-tailed_macaque\t9545\tMacaca nemestrina\tpig-tailed macaque\tGCF_043159975.1\tfalse",
                "common_marmoset\t9483\tCallithrix jacchus\tcommon marmoset\tGCF_049354715.1\tfalse",
                "Bolivian_squirrel_monkey\t27679\tSaimiri boliviensis\tBolivian squirrel monkey\tGCF_048565385.1\tfalse",
                "sunda_slow_loris\t9470\tNycticebus coucang\tsunda slow loris\tGCF_027406575.1\tfalse",
                "ring-tailed_lemur\t9447\tLemur catta\tring-tailed lemur\tGCF_020740605.2\tfalse",
                "Philippine_flying_lemur\t110931\tCynocephalus volans\tPhilippine flying lemur\tGCF_027409185.1\ttrue",
                "",
            ]
        )
    )
    log(f"Wrote manifest: {target}")


def sanitize_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    token = token.strip("._")
    return token or "species"


def collect_gcfs(args: argparse.Namespace) -> List[str]:
    gcfs: List[str] = []
    if args.inputfile:
        with Path(args.inputfile).open() as fh:
            for line in fh:
                text = line.strip()
                if not text or text.startswith("#"):
                    continue
                gcfs.append(text.split()[0])
    gcfs.extend(args.gcfs or [])
    seen = set()
    unique = []
    for gcf in gcfs:
        if gcf not in seen:
            seen.add(gcf)
            unique.append(gcf)
    if len(unique) < 2:
        raise ValueError("Provide at least two GCF accessions")
    return unique


def summary_for_gcf(gcf: str) -> dict:
    datasets = ROOT / "bin" / "datasets"
    require_file(datasets, "NCBI Datasets CLI")
    stdout = run_capture([
        str(datasets),
        "summary",
        "genome",
        "accession",
        gcf,
        "--as-json-lines",
    ])
    for line in stdout.splitlines():
        if line.strip():
            return json.loads(line)
    raise RuntimeError(f"No NCBI Datasets summary returned for {gcf}")


def organism_fields(record: dict) -> tuple[int, str]:
    organism = record.get("organism") or {}
    taxid = organism.get("taxId") or organism.get("tax_id") or organism.get("taxid")
    name = organism.get("organismName") or organism.get("organism_name") or organism.get("name") or ""
    if not taxid:
        raise ValueError(f"Missing organism taxid in record for {record.get('accession')}")
    return int(taxid), name


def unique_tokens(rows: List[dict]) -> None:
    counts: dict[str, int] = {}
    for row in rows:
        base = row["token"]
        counts[base] = counts.get(base, 0) + 1
        if counts[base] > 1:
            row["token"] = f"{base}_{row['taxid']}"
    counts.clear()
    for row in rows:
        base = row["token"]
        counts[base] = counts.get(base, 0) + 1
        if counts[base] > 1:
            suffix = re.sub(r"[^A-Za-z0-9]+", "_", row["gcf_accession"]).strip("_")
            row["token"] = f"{base}_{suffix}"


def cmd_manifest_from_gcf(args: argparse.Namespace) -> None:
    if args.download_tools:
        download_datasets_cli()
    gcfs = collect_gcfs(args)
    out = Path(args.output)
    if out.exists() and not args.force:
        raise FileExistsError(f"{out} already exists. Use --force to overwrite.")

    rows: List[dict] = []
    for gcf in gcfs:
        rec = summary_for_gcf(gcf)
        taxid, scientific_name = organism_fields(rec)
        if taxid in KNOWN_TAXID_LABELS:
            token, common_name = KNOWN_TAXID_LABELS[taxid]
        else:
            token = sanitize_token(scientific_name)
            common_name = scientific_name
        rows.append(
            {
                "token": token,
                "taxid": taxid,
                "scientific_name": scientific_name,
                "common_name": common_name,
                "gcf_accession": rec.get("accession") or gcf,
                "outgroup": "false",
            }
        )

    taxids = [r["taxid"] for r in rows]
    duplicates = sorted({t for t in taxids if taxids.count(t) > 1})
    if duplicates:
        raise ValueError(
            "Duplicate taxids are not supported because NCBI Gene orthology is taxid/GeneID-based: "
            f"{duplicates}"
        )
    unique_tokens(rows)

    if args.outgroup:
        matched = False
        for row in rows:
            if args.outgroup in {row["token"], row["gcf_accession"], str(row["taxid"])}:
                row["outgroup"] = "true"
                matched = True
        if not matched:
            raise ValueError(f"--outgroup did not match a token, taxid, or GCF accession: {args.outgroup}")

    if 9606 not in taxids:
        log("WARNING: manifest does not include taxid 9606. Human CDS matrices will be unavailable.")

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "token",
                "taxid",
                "scientific_name",
                "common_name",
                "gcf_accession",
                "outgroup",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    log(f"Wrote GCF-derived manifest with {len(rows)} assemblies: {out}")


def cmd_run(args: argparse.Namespace) -> None:
    require_file(RUN_DRIVER, "pipeline driver")
    if args.download_tools:
        download_datasets_cli()
    snapshot = ROOT / "reports" / "run_manifest_snapshot.tsv"
    if snapshot.exists() and DEFAULT_MANIFEST.exists() and not args.force:
        if snapshot.read_text() != DEFAULT_MANIFEST.read_text():
            raise RuntimeError(
                "config/species_manifest.tsv differs from the manifest used for existing outputs. "
                "Run with --force to rebuild all stages for the new GCF set."
            )
    cmd = [
        python_executable(),
        str(RUN_DRIVER),
        "--steps",
        args.steps,
        "--reference-taxid",
        str(args.reference_taxid),
    ]
    if args.force:
        cmd.append("--force")
    run(cmd)
    if DEFAULT_MANIFEST.exists():
        snapshot.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(DEFAULT_MANIFEST, snapshot)
    if args.with_matrices:
        if manifest_has_taxid(9606):
            matrix_args = argparse.Namespace(force=args.force)
            cmd_build_matrices(matrix_args)
        else:
            log("Skipping human CDS matrices because manifest does not include taxid 9606")


def cmd_build_matrices(args: argparse.Namespace) -> None:
    require_file(MATRIX_DRIVER, "human CDS matrix driver")
    cmd = [python_executable(), str(MATRIX_DRIVER)]
    if args.force:
        cmd.append("--force")
    run(cmd)


def cmd_align(args: argparse.Namespace) -> None:
    require_file(ALIGN_DRIVER, "MAFFT+PAL2NAL alignment driver")
    if args.mode != "mafft-pal2nal":
        raise ValueError(f"Unsupported alignment mode: {args.mode}")
    cmd = [
        python_executable(),
        str(ALIGN_DRIVER),
        "--input-dir",
        args.input_dir,
        "--output-dir",
        args.output_dir,
        "--mafft",
        args.mafft,
        "--pal2nal",
        args.pal2nal,
        "--threads-per-mafft",
        str(args.threads_per_mafft),
        "--jobs",
        str(args.jobs),
        "--codon-table",
        args.codon_table,
    ]
    if args.limit is not None:
        cmd.extend(["--limit", str(args.limit)])
    if args.symbols:
        cmd.append("--symbols")
        cmd.extend(args.symbols)
    for token in args.map_token or []:
        cmd.extend(["--map-token", token])
    if args.force:
        cmd.append("--force")
    run(cmd)


def read_tsv(path: Path) -> List[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def manifest_has_taxid(taxid: int) -> bool:
    if not DEFAULT_MANIFEST.exists():
        return False
    for row in read_tsv(DEFAULT_MANIFEST):
        try:
            if int(row["taxid"]) == taxid:
                return True
        except Exception:
            continue
    return False


def count_fasta_records(path: Path) -> tuple[int, set[str]]:
    count = 0
    headers: set[str] = set()
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
                headers.add(line[1:].strip().split()[0])
    return count, headers


def count_gzip_tsv_rows(path: Path) -> int:
    with gzip.open(path, "rt") as fh:
        next(fh, None)
        return sum(1 for _ in fh)


def verify_python_sources() -> None:
    for path in [RUN_DRIVER, MATRIX_DRIVER, ALIGN_DRIVER, ROOT / "refseq2cds.py"]:
        require_file(path, "Python source")
        run([python_executable(), "-m", "py_compile", str(path)])


def manifest_tokens() -> set[str]:
    manifest = DEFAULT_MANIFEST
    require_file(manifest, "species manifest")
    rows = read_tsv(manifest)
    return {row["token"] for row in rows}


def verify_fastas(full: bool) -> dict:
    manifest_path = ROOT / "fastas" / "manifest.tsv"
    require_file(manifest_path, "FASTA manifest")
    rows = read_tsv(manifest_path)
    fasta_files = sorted((ROOT / "fastas").glob("*.fasta"))
    expected_tokens = manifest_tokens()
    expected_count = len(expected_tokens)
    if len(rows) != len(fasta_files):
        raise AssertionError(f"FASTA manifest rows {len(rows)} != FASTA files {len(fasta_files)}")
    check_rows = rows if full else rows[:200]
    bad = []
    for row in check_rows:
        fasta_path = ROOT / row["fasta_path"]
        count, headers = count_fasta_records(fasta_path)
        if count != expected_count or headers != expected_tokens:
            bad.append((row["fasta_path"], count, sorted(headers ^ expected_tokens)))
    if bad:
        raise AssertionError(f"Bad FASTA files: {bad[:5]}")
    return {
        "manifest_rows": len(rows),
        "fasta_files": len(fasta_files),
        "checked": len(check_rows),
        "expected_sequence_count": expected_count,
    }


def verify_matrices(mode: str) -> dict:
    manifest_path = ROOT / "human_cds_matrices" / "manifest.tsv"
    require_file(manifest_path, "human CDS matrix manifest")
    rows = read_tsv(manifest_path)
    matrix_files = sorted((ROOT / "human_cds_matrices").glob("*.human_cds_genomic_matrix.tsv.gz"))
    failed_path = ROOT / "human_cds_matrices" / "failed.tsv"
    require_file(failed_path, "human CDS matrix failure log")
    failed_rows = max(0, sum(1 for _ in failed_path.open()) - 1)
    if failed_rows != 0:
        raise AssertionError(f"Matrix failures are present: {failed_rows}")
    if len(rows) != len(matrix_files):
        raise AssertionError(f"Matrix manifest rows {len(rows)} != matrix files {len(matrix_files)}")
    if mode == "none":
        checked = 0
    elif mode == "sample":
        interesting = {"BRCA1", "TP53", "A1BG", "APOB"}
        sample = [r for r in rows if r.get("human_symbol") in interesting]
        if not sample:
            sample = rows[:20]
        checked = len(sample)
        for row in sample:
            actual = count_gzip_tsv_rows(ROOT / row["matrix_path"])
            expected = int(row["matrix_rows"])
            if actual != expected:
                raise AssertionError(f"{row['matrix_path']} rows {actual} != manifest {expected}")
    elif mode == "full":
        checked = len(rows)
        for row in rows:
            actual = count_gzip_tsv_rows(ROOT / row["matrix_path"])
            expected = int(row["matrix_rows"])
            if actual != expected:
                raise AssertionError(f"{row['matrix_path']} rows {actual} != manifest {expected}")
    else:
        raise ValueError(f"Unknown matrix check mode: {mode}")
    return {
        "manifest_rows": len(rows),
        "matrix_files": len(matrix_files),
        "failed_rows": failed_rows,
        "row_count_checked": checked,
    }


def verify_required_reports() -> dict:
    summary = ROOT / "reports" / "summary.json"
    matrix_summary = ROOT / "reports" / "human_cds_position_matrices.summary.json"
    require_file(summary, "pipeline summary")
    require_file(matrix_summary, "matrix summary")
    return {
        "summary": json.loads(summary.read_text()),
        "matrix_summary": json.loads(matrix_summary.read_text()),
    }


def cmd_verify(args: argparse.Namespace) -> None:
    verify_python_sources()
    fasta_stats = verify_fastas(full=args.full)
    matrix_stats = verify_matrices(mode=args.matrix_rows)
    report_stats = verify_required_reports()
    result = {
        "status": "pass",
        "fastas": fasta_stats,
        "matrices": matrix_stats,
        "summary_fastas": report_stats["summary"].get("fastas", {}),
        "summary_matrix_files": report_stats["matrix_summary"].get("matrix_files"),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


def cmd_summary(_args: argparse.Namespace) -> None:
    verify_required_reports()
    for path in [
        ROOT / "reports" / "summary.json",
        ROOT / "reports" / "human_cds_position_matrices.summary.json",
    ]:
        print(f"\n# {path.relative_to(ROOT)}")
        print(path.read_text())


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="refseq2cds",
        description="Assembly-exact NCBI RefSeq singleton ortholog CDS FASTA and human CDS-position matrix builder.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("init-manifest", help="Write the default 14-species manifest")
    p.add_argument("-o", "--output", default=str(DEFAULT_MANIFEST))
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_init_manifest)

    p = sub.add_parser("manifest-from-gcf", help="Create config/species_manifest.tsv from GCF accessions")
    p.add_argument("gcfs", nargs="*", help="GCF accessions")
    p.add_argument("-i", "--inputfile", help="File with one GCF accession per line")
    p.add_argument("-o", "--output", default=str(DEFAULT_MANIFEST))
    p.add_argument("--outgroup", help="Optional outgroup token, taxid, or GCF accession to mark in the manifest")
    p.add_argument("--download-tools", action="store_true", help="Download NCBI Datasets CLI before querying")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_manifest_from_gcf)

    p = sub.add_parser("download-tools", help="Download NCBI Datasets CLI binaries into ./bin")
    p.set_defaults(func=lambda _args: download_datasets_cli())

    p = sub.add_parser("run", help="Run the CDS pipeline")
    p.add_argument("--steps", default="all", help="Pipeline steps, comma-separated, or all")
    p.add_argument(
        "--reference-taxid",
        type=int,
        default=9606,
        help="TaxID whose gene symbol is used for family IDs and FASTA filenames (default: 9606)",
    )
    p.add_argument("--force", action="store_true")
    p.add_argument("--download-tools", action="store_true")
    p.add_argument("--with-matrices", action="store_true", help="Build human CDS genomic matrices after FASTA generation")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("build-matrices", help="Build human CDS genomic-position matrices")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_build_matrices)

    p = sub.add_parser("align", help="Build tree-free codon alignments")
    p.add_argument(
        "--mode",
        choices=["mafft-pal2nal"],
        default="mafft-pal2nal",
        help="Alignment backend (default: mafft-pal2nal)",
    )
    p.add_argument("--input-dir", default=str(ROOT / "fastas"), help="Directory containing refseq2cds FASTA files")
    p.add_argument(
        "--output-dir",
        default=str(ROOT / "alignments" / "mafft_pal2nal"),
        help="Alignment output directory",
    )
    p.add_argument("--mafft", default="mafft", help="MAFFT executable path or command name")
    p.add_argument("--pal2nal", default="pal2nal.pl", help="PAL2NAL executable path or command name")
    p.add_argument("--threads-per-mafft", type=int, default=1, help="Threads passed to each MAFFT process")
    p.add_argument("--jobs", type=int, default=1, help="Number of families to align concurrently")
    p.add_argument("--limit", type=int, help="Only process the first N selected FASTA files")
    p.add_argument("--symbols", nargs="*", help="Only process selected symbols/stems; comma-separated values are accepted")
    p.add_argument("--map-token", action="append", default=[], help="Write alignment-to-CDS codon map for this token")
    p.add_argument("--codon-table", choices=["universal", "vmitochondria"], default="universal")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_align)

    p = sub.add_parser("verify", help="Verify code and generated outputs")
    p.add_argument("--full", action="store_true", help="Check every FASTA instead of a 200-file sample")
    p.add_argument(
        "--matrix-rows",
        choices=["none", "sample", "full"],
        default="sample",
        help="How deeply to count gzip matrix rows",
    )
    p.set_defaults(func=cmd_verify)

    p = sub.add_parser("summary", help="Print generated summary JSON files")
    p.set_defaults(func=cmd_summary)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
