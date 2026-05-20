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
import tempfile
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
    candidates = [
        Path.cwd().resolve(),
        MODULE_ROOT,
        Path(sys.prefix).resolve() / "share" / "refseq2cds",
        Path(sys.exec_prefix).resolve() / "share" / "refseq2cds",
    ]
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
VARIANT_DRIVER = ROOT / "workflow" / "scripts" / "call_target_specific_variants.py"
VERSION = "0.1.5"

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

GCF_ACCESSION_RE = re.compile(r"(GCF_\d+\.\d+)")


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


def default_manifest_path() -> Path:
    local = Path.cwd() / "config" / "species_manifest.tsv"
    if local.exists():
        return local.resolve()
    return DEFAULT_MANIFEST


def resolve_tool(name: str) -> str:
    """Resolve an executable from explicit path, ./bin, package bin, or PATH."""
    path = Path(name)
    if path.exists():
        return str(path.resolve())
    if path.parent != Path("."):
        return str(path)
    for candidate in [Path.cwd() / "bin" / name, ROOT / "bin" / name]:
        if candidate.exists():
            return str(candidate.resolve())
    found = shutil.which(name)
    if found:
        return found
    return name


def download_datasets_cli() -> None:
    bin_dir = Path.cwd() / "bin"
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


def normalize_gcf_accession(value: str) -> str:
    """Return the base RefSeq assembly accession from common user inputs.

    Some planning tables append RefSeq annotation-release labels such as
    ``-RS_2025_03`` or include a trailing slash. The NCBI Datasets genome
    accession command expects the base ``GCF_<digits>.<version>`` accession.
    """
    match = GCF_ACCESSION_RE.search(value.strip())
    if not match:
        raise ValueError(f"Could not parse a RefSeq GCF accession from: {value!r}")
    return match.group(1)


def collect_gcfs(args: argparse.Namespace) -> List[str]:
    gcfs: List[str] = []
    if args.inputfile:
        with Path(args.inputfile).open() as fh:
            for line in fh:
                text = line.strip()
                if not text or text.startswith("#"):
                    continue
                gcfs.append(normalize_gcf_accession(text.split()[0]))
    gcfs.extend(normalize_gcf_accession(gcf) for gcf in (args.gcfs or []))
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
    datasets = resolve_tool("datasets")
    if not shutil.which(datasets) and not Path(datasets).exists():
        raise FileNotFoundError(
            "Missing NCBI Datasets CLI. Run `refseq2cds download-tools` or install ncbi-datasets-cli."
        )
    stdout = run_capture([
        datasets,
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
    datasets = resolve_tool(args.datasets)
    manifest = Path(args.manifest).resolve() if args.manifest else default_manifest_path()
    output_root = Path(args.output_root).resolve() if args.output_root else Path.cwd().resolve()
    snapshot = output_root / "reports" / "run_manifest_snapshot.tsv"
    config_snapshot = output_root / "reports" / "run_config_snapshot.json"
    if snapshot.exists() and manifest.exists() and not args.force:
        if snapshot.read_text() != manifest.read_text():
            raise RuntimeError(
                f"{manifest} differs from the manifest used for existing outputs. "
                "Run with --force to rebuild all stages for the new GCF set."
            )
    current_config = {
        "orthology_mode": args.orthology_mode,
        "reference_taxid": int(args.reference_taxid),
        "reference_symbol": args.reference_symbol or "",
        "reference_gene_id": args.reference_gene_id,
        "min_sequences": int(args.min_sequences) if args.min_sequences is not None else None,
        "exclude_reference": bool(args.exclude_reference),
    }
    if config_snapshot.exists() and not args.force:
        previous_config = json.loads(config_snapshot.read_text())
        if previous_config != current_config:
            raise RuntimeError(
                f"{config_snapshot} differs from the requested run configuration. "
                "Run with --force to rebuild outputs for the new orthology/reference settings."
            )
    cmd = [
        python_executable(),
        str(RUN_DRIVER),
        "--steps",
        args.steps,
        "--orthology-mode",
        args.orthology_mode,
        "--reference-taxid",
        str(args.reference_taxid),
        "--manifest",
        str(manifest),
        "--output-root",
        str(output_root),
        "--datasets",
        datasets,
    ]
    if args.input_root:
        cmd.extend(["--input-root", args.input_root])
    if args.offline:
        cmd.append("--offline")
    if args.force:
        cmd.append("--force")
    if args.reference_symbol:
        cmd.extend(["--reference-symbol", args.reference_symbol])
    if args.reference_gene_id is not None:
        cmd.extend(["--reference-gene-id", str(args.reference_gene_id)])
    if args.min_sequences is not None:
        cmd.extend(["--min-sequences", str(args.min_sequences)])
    if args.exclude_reference:
        cmd.append("--exclude-reference")
    run(cmd)
    if manifest.exists():
        snapshot.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(manifest, snapshot)
    config_snapshot.parent.mkdir(parents=True, exist_ok=True)
    config_snapshot.write_text(json.dumps(current_config, indent=2, sort_keys=True) + "\n")
    if args.with_matrices:
        if manifest_has_taxid(9606, manifest):
            matrix_args = argparse.Namespace(
                force=args.force,
                manifest=str(manifest),
                input_root=args.input_root,
                output_root=str(output_root),
            )
            cmd_build_matrices(matrix_args)
        else:
            log("Skipping human CDS matrices because manifest does not include taxid 9606")


def cmd_build_matrices(args: argparse.Namespace) -> None:
    require_file(MATRIX_DRIVER, "human CDS matrix driver")
    output_root = Path(args.output_root).resolve() if args.output_root else Path.cwd().resolve()
    manifest = Path(args.manifest).resolve() if args.manifest else default_manifest_path()
    cmd = [
        python_executable(),
        str(MATRIX_DRIVER),
        "--manifest",
        str(manifest),
    ]
    if args.input_root:
        cmd.extend(["--input-root", args.input_root])
    cmd.extend(["--output-root", str(output_root)])
    if args.force:
        cmd.append("--force")
    run(cmd)


def cmd_align(args: argparse.Namespace) -> None:
    require_file(ALIGN_DRIVER, "MAFFT+PAL2NAL alignment driver")
    if args.mode != "mafft-pal2nal":
        raise ValueError(f"Unsupported alignment mode: {args.mode}")
    input_dir = Path(args.input_dir).resolve()
    fasta_manifest = (
        Path(args.fasta_manifest).resolve()
        if args.fasta_manifest
        else (input_dir / "manifest.tsv").resolve()
    )
    default_species_manifest = input_dir.parent / "reports" / "run_manifest_snapshot.tsv"
    species_manifest = (
        Path(args.species_manifest).resolve()
        if args.species_manifest
        else (default_species_manifest.resolve() if default_species_manifest.exists() else default_manifest_path())
    )
    cmd = [
        python_executable(),
        str(ALIGN_DRIVER),
        "--input-dir",
        str(input_dir),
        "--output-dir",
        args.output_dir,
        "--fasta-manifest",
        str(fasta_manifest),
        "--species-manifest",
        str(species_manifest),
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


def cmd_variants(args: argparse.Namespace) -> None:
    require_file(VARIANT_DRIVER, "target-specific variant driver")
    cmd = [
        python_executable(),
        str(VARIANT_DRIVER),
        "--alignment-dir",
        args.alignment_dir,
        "--matrix-dir",
        args.matrix_dir,
        "--output-dir",
        args.output_dir,
        "--coordinate-reference-token",
        args.coordinate_reference_token,
        "--target-state-mode",
        args.target_state_mode,
        "--min-target-non-gap",
        args.min_target_non_gap,
        "--min-target-informative",
        args.min_target_informative if args.min_target_informative is not None else args.min_target_non_gap,
        "--max-target-gap-fraction-for-substitution",
        str(args.max_target_gap_fraction_for_substitution),
        "--min-background-non-gap",
        str(args.min_background_non_gap),
        "--min-background-informative",
        args.min_background_informative if args.min_background_informative is not None else str(args.min_background_non_gap),
        "--min-background-gap-fraction-for-target-non-gap-event",
        str(args.min_background_gap_fraction_for_target_non_gap_event),
        "--min-background-non-gap-fraction-for-target-gap-event",
        str(args.min_background_non_gap_fraction_for_target_gap_event),
        "--min-target-non-gap-fraction",
        str(args.min_target_non_gap_fraction),
        "--min-target-gap-fraction",
        str(args.min_target_gap_fraction),
        "--bed-mode",
        args.bed_mode,
        "--codon-table",
        args.codon_table,
    ]
    if args.codon_map_dir:
        cmd.extend(["--codon-map-dir", args.codon_map_dir])
    if args.reference_token:
        cmd.extend(["--reference-token", args.reference_token])
    for token in args.target_token or []:
        cmd.extend(["--target-token", token])
    if args.target_tokens_file:
        cmd.extend(["--target-tokens-file", args.target_tokens_file])
    for token in args.outgroup_token or []:
        cmd.extend(["--outgroup-token", token])
    if args.outgroup_tokens_file:
        cmd.extend(["--outgroup-tokens-file", args.outgroup_tokens_file])
    for token in args.exclude_token or []:
        cmd.extend(["--exclude-token", token])
    if args.exclude_tokens_file:
        cmd.extend(["--exclude-tokens-file", args.exclude_tokens_file])
    if args.symbols:
        cmd.append("--symbols")
        cmd.extend(args.symbols)
    if args.symbols_file:
        cmd.extend(["--symbols-file", args.symbols_file])
    if args.limit is not None:
        cmd.extend(["--limit", str(args.limit)])
    if args.force:
        cmd.append("--force")
    run(cmd)


def read_tsv(path: Path) -> List[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def manifest_has_taxid(taxid: int, manifest: Path = DEFAULT_MANIFEST) -> bool:
    if not manifest.exists():
        return False
    for row in read_tsv(manifest):
        try:
            if int(row["taxid"]) == taxid:
                return True
        except Exception:
            continue
    return False


def cmd_test(args: argparse.Namespace) -> None:
    if args.example != "mini":
        raise ValueError(f"Unsupported example: {args.example}")
    example_dir = ROOT / "examples" / "mini"
    manifest = example_dir / "manifest.tsv"
    input_root = example_dir / "inputs"
    require_file(manifest, "mini manifest")
    require_file(input_root / "ncbi_bulk" / "gene_orthologs.mini.tsv", "mini gene_orthologs fixture")
    if args.output_root:
        output_root = Path(args.output_root).resolve()
        if output_root.exists() and args.force:
            shutil.rmtree(output_root)
        output_root.mkdir(parents=True, exist_ok=True)
    else:
        output_root = Path(tempfile.mkdtemp(prefix="refseq2cds_mini_"))
    run_args = argparse.Namespace(
        steps="all",
        orthology_mode="strict_singleton",
        reference_taxid=9606,
        reference_symbol=None,
        reference_gene_id=None,
        min_sequences=2,
        exclude_reference=False,
        manifest=str(manifest),
        input_root=str(input_root),
        output_root=str(output_root),
        datasets=resolve_tool("datasets"),
        offline=True,
        force=True,
        download_tools=False,
        with_matrices=True,
    )
    cmd_run(run_args)
    summary = json.loads((output_root / "reports" / "summary.json").read_text())
    matrix_summary = json.loads((output_root / "reports" / "human_cds_position_matrices.summary.json").read_text())
    fasta_count = int(summary.get("fastas", {}).get("count", -1))
    matrix_count = int(matrix_summary.get("matrix_files", -1))
    matrix_failed = int(matrix_summary.get("failed", -1))
    rejected = summary.get("rejected_reasons", {})
    if fasta_count != 2:
        raise AssertionError(f"mini expected 2 FASTA files, observed {fasta_count}")
    if matrix_count != 2 or matrix_failed != 0:
        raise AssertionError(f"mini expected 2 successful matrices, observed count={matrix_count}, failed={matrix_failed}")
    if int(rejected.get("component_size_not_expected_species_count", 0)) < 2:
        raise AssertionError(f"mini expected missing/paralog rejection counts, observed {rejected}")
    result = {
        "status": "pass",
        "example": "mini",
        "output_root": str(output_root),
        "fastas": fasta_count,
        "human_matrices": matrix_count,
        "rejected_reasons": rejected,
    }
    print(json.dumps(result, indent=2, sort_keys=True))


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
    for path in [RUN_DRIVER, MATRIX_DRIVER, ALIGN_DRIVER, VARIANT_DRIVER, Path(__file__).resolve()]:
        require_file(path, "Python source")
        run([python_executable(), "-m", "py_compile", str(path)])


def manifest_tokens(manifest: Path) -> set[str]:
    require_file(manifest, "species manifest")
    rows = read_tsv(manifest)
    return {row["token"] for row in rows}


def verify_fastas(full: bool, output_root: Path, species_manifest: Path) -> dict:
    manifest_path = output_root / "fastas" / "manifest.tsv"
    require_file(manifest_path, "FASTA manifest")
    rows = read_tsv(manifest_path)
    fasta_files = sorted(
        p
        for p in (output_root / "fastas").glob("*.fasta")
        if not p.name.startswith("._")
    )
    expected_tokens = manifest_tokens(species_manifest)
    expected_count = len(expected_tokens)
    if len(rows) != len(fasta_files):
        raise AssertionError(f"FASTA manifest rows {len(rows)} != FASTA files {len(fasta_files)}")
    check_rows = rows if full else rows[:200]
    bad = []
    for row in check_rows:
        fasta_path = output_root / row["fasta_path"]
        count, headers = count_fasta_records(fasta_path)
        mode = row.get("orthology_mode", "strict_singleton")
        if mode == "reference_gene_1to1_present_species":
            manifest_count = int(row.get("sequence_count", count))
            if count != manifest_count or not headers.issubset(expected_tokens):
                bad.append((row["fasta_path"], count, sorted(headers - expected_tokens)))
        elif count != expected_count or headers != expected_tokens:
            bad.append((row["fasta_path"], count, sorted(headers ^ expected_tokens)))
    if bad:
        raise AssertionError(f"Bad FASTA files: {bad[:5]}")
    return {
        "manifest_rows": len(rows),
        "fasta_files": len(fasta_files),
        "checked": len(check_rows),
        "expected_sequence_count": expected_count,
    }


def verify_matrices(mode: str, output_root: Path) -> dict:
    """Verify human CDS-to-genome matrices.

    Matrix generation is optional. Passing ``--matrix-rows none`` means the
    caller wants a FASTA/report-only verification, so matrix files are not
    required in that mode.
    """
    if mode == "none":
        return {
            "status": "skipped",
            "manifest_rows": 0,
            "matrix_files": 0,
            "failed_rows": 0,
            "row_count_checked": 0,
        }
    manifest_path = output_root / "human_cds_matrices" / "manifest.tsv"
    require_file(manifest_path, "human CDS matrix manifest")
    rows = read_tsv(manifest_path)
    matrix_files = sorted(
        p
        for p in (output_root / "human_cds_matrices").glob("*.human_cds_genomic_matrix.tsv.gz")
        if not p.name.startswith("._")
    )
    failed_path = output_root / "human_cds_matrices" / "failed.tsv"
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
            actual = count_gzip_tsv_rows(output_root / row["matrix_path"])
            expected = int(row["matrix_rows"])
            if actual != expected:
                raise AssertionError(f"{row['matrix_path']} rows {actual} != manifest {expected}")
    elif mode == "full":
        checked = len(rows)
        for row in rows:
            actual = count_gzip_tsv_rows(output_root / row["matrix_path"])
            expected = int(row["matrix_rows"])
            if actual != expected:
                raise AssertionError(f"{row['matrix_path']} rows {actual} != manifest {expected}")
    else:
        raise ValueError(f"Unknown matrix check mode: {mode}")
    return {
        "status": "pass",
        "manifest_rows": len(rows),
        "matrix_files": len(matrix_files),
        "failed_rows": failed_rows,
        "row_count_checked": checked,
    }


def verify_required_reports(output_root: Path, require_matrices: bool = True) -> dict:
    summary = output_root / "reports" / "summary.json"
    matrix_summary = output_root / "reports" / "human_cds_position_matrices.summary.json"
    require_file(summary, "pipeline summary")
    result = {"summary": json.loads(summary.read_text())}
    if require_matrices:
        require_file(matrix_summary, "matrix summary")
        result["matrix_summary"] = json.loads(matrix_summary.read_text())
    else:
        result["matrix_summary"] = None
    return result


def count_data_rows(path: Path) -> int:
    require_file(path, "TSV file")
    with path.open() as fh:
        return max(0, sum(1 for _ in fh) - 1)


def count_bed_rows(path: Path) -> int:
    require_file(path, "BED file")
    with path.open() as fh:
        return sum(1 for line in fh if line.strip())


def manifest_bed_paths(variant_dir: Path) -> List[Path]:
    manifest = variant_dir / "manifest.tsv"
    paths: List[Path] = []
    with manifest.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            for path_text in row.get("bed_paths", "").split(";"):
                if path_text:
                    paths.append((variant_dir / path_text).resolve())
    return paths


def verify_alignments(alignment_dir: Optional[Path]) -> dict:
    """Verify optional MAFFT+PAL2NAL outputs when requested.

    The align stage consumes ``fastas/`` and writes a self-contained directory
    with ``codon/``, ``maps/``, ``manifest.tsv``, ``failed.tsv``, and
    ``summary.json``. This check deliberately stays lightweight: it verifies
    counts and failure status without rereading every alignment sequence.
    """
    if alignment_dir is None:
        return {"status": "skipped"}
    require_file(alignment_dir / "summary.json", "alignment summary")
    require_file(alignment_dir / "manifest.tsv", "alignment manifest")
    require_file(alignment_dir / "failed.tsv", "alignment failure log")
    summary = json.loads((alignment_dir / "summary.json").read_text())
    manifest_rows = count_data_rows(alignment_dir / "manifest.tsv")
    failed_rows = count_data_rows(alignment_dir / "failed.tsv")
    codon_dir = alignment_dir / "codon"
    codon_files = sorted(p for p in codon_dir.glob("*.codon.fasta") if not p.name.startswith("._"))
    if failed_rows != 0 or int(summary.get("failed", 0)) != 0:
        raise AssertionError(f"Alignment failures are present: failed.tsv={failed_rows}, summary={summary.get('failed')}")
    if manifest_rows != len(codon_files):
        raise AssertionError(f"Alignment manifest rows {manifest_rows} != codon alignment files {len(codon_files)}")
    if int(summary.get("passed", manifest_rows)) != manifest_rows:
        raise AssertionError(f"Alignment summary passed {summary.get('passed')} != manifest rows {manifest_rows}")
    return {
        "status": "pass",
        "manifest_rows": manifest_rows,
        "codon_alignment_files": len(codon_files),
        "failed_rows": failed_rows,
        "summary_passed": int(summary.get("passed", manifest_rows)),
    }


def verify_variants(variant_dir: Optional[Path]) -> dict:
    """Verify optional variant-call outputs when requested."""
    if variant_dir is None:
        return {"status": "skipped"}
    require_file(variant_dir / "summary.json", "variant summary")
    require_file(variant_dir / "manifest.tsv", "variant manifest")
    summary = json.loads((variant_dir / "summary.json").read_text())
    manifest_rows = count_data_rows(variant_dir / "manifest.tsv")
    if summary.get("status") != "pass" or int(summary.get("families_failed", 0)) != 0:
        raise AssertionError(
            f"Variant failures are present: status={summary.get('status')}, "
            f"families_failed={summary.get('families_failed')}"
        )
    if int(summary.get("families_processed", manifest_rows)) != manifest_rows:
        raise AssertionError(f"Variant summary families_processed {summary.get('families_processed')} != manifest rows {manifest_rows}")
    bed_paths = manifest_bed_paths(variant_dir)
    if not bed_paths and (variant_dir / "bed").exists():
        bed_paths = [
            path.resolve()
            for path in (variant_dir / "bed").glob("*.bed")
            if not path.name.startswith("._") and path.name != "merged.bed"
        ]
    bed_rows = sum(count_bed_rows(path) for path in bed_paths)
    if int(summary.get("bed_rows", bed_rows)) != bed_rows:
        raise AssertionError(f"Variant summary bed_rows {summary.get('bed_rows')} != per-gene BED rows {bed_rows}")
    merged_bed_path = summary.get("merged_bed_path")
    if merged_bed_path:
        merged_path = (variant_dir / merged_bed_path).resolve()
        merged_rows = count_bed_rows(merged_path)
        if int(summary.get("merged_bed_rows", merged_rows)) != merged_rows:
            raise AssertionError(
                f"Variant summary merged_bed_rows {summary.get('merged_bed_rows')} != merged BED rows {merged_rows}"
            )
        if merged_rows != bed_rows:
            raise AssertionError(f"Merged BED rows {merged_rows} != per-gene BED rows {bed_rows}")
    else:
        merged_rows = 0
    return {
        "status": "pass",
        "manifest_rows": manifest_rows,
        "families_processed": int(summary.get("families_processed", manifest_rows)),
        "bed_rows": bed_rows,
        "merged_bed_rows": merged_rows,
    }


def resolve_verify_subdir(value: Optional[str], output_root: Path) -> Optional[Path]:
    if not value:
        return None
    path = Path(value)
    return path.resolve() if path.is_absolute() else (output_root / path).resolve()


def cmd_verify(args: argparse.Namespace) -> None:
    output_root = Path(args.output_root).resolve() if args.output_root else Path.cwd().resolve()
    species_manifest = Path(args.manifest).resolve() if args.manifest else default_manifest_path()
    verify_python_sources()
    fasta_stats = verify_fastas(full=args.full, output_root=output_root, species_manifest=species_manifest)
    matrix_stats = verify_matrices(mode=args.matrix_rows, output_root=output_root)
    report_stats = verify_required_reports(output_root, require_matrices=args.matrix_rows != "none")
    alignment_dir = resolve_verify_subdir(args.alignment_dir, output_root)
    variant_dir = resolve_verify_subdir(args.variant_dir, output_root)
    alignment_stats = verify_alignments(alignment_dir)
    variant_stats = verify_variants(variant_dir)
    result = {
        "status": "pass",
        "fastas": fasta_stats,
        "matrices": matrix_stats,
        "alignments": alignment_stats,
        "variants": variant_stats,
        "summary_fastas": report_stats["summary"].get("fastas", {}),
        "summary_matrix_files": (
            report_stats["matrix_summary"].get("matrix_files")
            if report_stats["matrix_summary"]
            else None
        ),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


def cmd_summary(args: argparse.Namespace) -> None:
    output_root = Path(args.output_root).resolve() if args.output_root else Path.cwd().resolve()
    verify_required_reports(output_root)
    for path in [
        output_root / "reports" / "summary.json",
        output_root / "reports" / "human_cds_position_matrices.summary.json",
    ]:
        print(f"\n# {path.relative_to(output_root)}")
        print(path.read_text())


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="refseq2cds",
        description="Assembly-exact RefSeq ortholog CDS FASTA, coordinate matrix, codon alignment, and variant workflow.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("init-manifest", help="Write the default 14-species manifest")
    p.add_argument("-o", "--output", default=str(Path("config") / "species_manifest.tsv"))
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_init_manifest)

    p = sub.add_parser("manifest-from-gcf", help="Create config/species_manifest.tsv from GCF accessions")
    p.add_argument("gcfs", nargs="*", help="GCF accessions")
    p.add_argument("-i", "--inputfile", help="File with one GCF accession per line")
    p.add_argument("-o", "--output", default=str(Path("config") / "species_manifest.tsv"))
    p.add_argument("--outgroup", help="Optional outgroup token, taxid, or GCF accession to mark in the manifest")
    p.add_argument("--download-tools", action="store_true", help="Download NCBI Datasets CLI before querying")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_manifest_from_gcf)

    p = sub.add_parser("download-tools", help="Download NCBI Datasets CLI binaries into ./bin")
    p.set_defaults(func=lambda _args: download_datasets_cli())

    p = sub.add_parser("run", help="Build assembly-exact ortholog CDS FASTA files and optional human matrices")
    p.add_argument("--steps", default="all", help="Pipeline stages to run, comma-separated, or all")
    p.add_argument(
        "--orthology-mode",
        choices=["strict_singleton", "reference_gene_1to1_present_species"],
        default="strict_singleton",
        help="Orthology mode: strict all-species singleton extraction or reference-gene present-species 1:1 extraction",
    )
    p.add_argument(
        "--reference-taxid",
        type=int,
        default=9606,
        help="TaxID whose gene symbol is used for family IDs and FASTA filenames (default: 9606)",
    )
    p.add_argument("--reference-symbol", help="Reference-species gene symbol for reference_gene_1to1_present_species mode")
    p.add_argument("--reference-gene-id", type=int, help="Reference-species NCBI GeneID for reference_gene_1to1_present_species mode")
    p.add_argument("--min-sequences", type=int, default=2, help="Minimum sequences required in reference-gene mode")
    p.add_argument("--exclude-reference", action="store_true", help="Exclude the reference sequence in reference-gene mode")
    p.add_argument(
        "--manifest",
        help="Species manifest TSV; defaults to ./config/species_manifest.tsv if present, otherwise the packaged default manifest",
    )
    p.add_argument("--input-root", help="Input root containing reusable ncbi_bulk/ and assembly_packages/")
    p.add_argument("--output-root", help="Run output root for indexes, FASTA files, reports, and matrices")
    p.add_argument("--datasets", default="datasets", help="NCBI Datasets CLI executable")
    p.add_argument("--offline", action="store_true", help="Use local input-root fixtures; do not download NCBI data")
    p.add_argument("--force", action="store_true")
    p.add_argument("--download-tools", action="store_true")
    p.add_argument("--with-matrices", action="store_true", help="Build human CDS genomic matrices after FASTA generation")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("build-matrices", help="Build human CDS genomic-position matrices")
    p.add_argument(
        "--manifest",
        help="Species manifest TSV; defaults to ./config/species_manifest.tsv if present, otherwise the packaged default manifest",
    )
    p.add_argument("--input-root", help="Input root containing assembly_packages/")
    p.add_argument("--output-root", help="Output root containing selection/ and fastas/")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_build_matrices)

    p = sub.add_parser("align", help="Build MAFFT+PAL2NAL codon alignments from fastas/")
    p.add_argument(
        "--mode",
        choices=["mafft-pal2nal"],
        default="mafft-pal2nal",
        help="Alignment backend (default: mafft-pal2nal)",
    )
    p.add_argument("--input-dir", default=str(Path.cwd() / "fastas"), help="Input FASTA directory from `refseq2cds run`")
    p.add_argument(
        "--output-dir",
        default=str(Path.cwd() / "alignments" / "mafft_pal2nal"),
        help="Alignment output directory consumed by `refseq2cds variants`",
    )
    p.add_argument("--mafft", default="mafft", help="MAFFT executable path or command name")
    p.add_argument("--pal2nal", default="pal2nal.pl", help="PAL2NAL executable path or command name")
    p.add_argument(
        "--fasta-manifest",
        help="FASTA manifest TSV; defaults to <input-dir>/manifest.tsv",
    )
    p.add_argument(
        "--species-manifest",
        help="Species manifest TSV; defaults to <input-dir>/../reports/run_manifest_snapshot.tsv when present",
    )
    p.add_argument("--threads-per-mafft", type=int, default=1, help="Threads passed to each MAFFT process")
    p.add_argument("--jobs", type=int, default=1, help="Number of families to align concurrently")
    p.add_argument("--limit", type=int, help="Only process the first N selected FASTA files")
    p.add_argument("--symbols", nargs="*", help="Only process selected symbols/stems; comma-separated values are accepted")
    p.add_argument("--map-token", action="append", default=[], help="Write alignment-to-CDS codon map for this token; required for later BED mapping")
    p.add_argument("--codon-table", choices=["universal", "vmitochondria"], default="universal")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_align)

    p = sub.add_parser("variants", help="Call target-specific amino acid variants and map coordinateable variants to BED")
    p.add_argument(
        "--alignment-dir",
        default=str(Path.cwd() / "alignments" / "mafft_pal2nal"),
        help="Directory containing codon-aware alignments, or an aligner output directory with a codon/ subdirectory",
    )
    p.add_argument(
        "--codon-map-dir",
        help="Directory containing {SYMBOL}.{TOKEN}.codon_map.tsv.gz files; defaults to maps near --alignment-dir",
    )
    p.add_argument("--matrix-dir", default=str(Path.cwd() / "human_cds_matrices"))
    p.add_argument("--output-dir", default=str(Path.cwd() / "variants"))
    p.add_argument(
        "--coordinate-reference-token",
        default="human",
        help="Token used only for CDS/genome coordinate mapping, not as a privileged event comparator",
    )
    p.add_argument("--reference-token", help="Deprecated alias for --coordinate-reference-token")
    p.add_argument("--target-token", action="append", default=[], help="Target token; may be used more than once")
    p.add_argument("--target-tokens-file", help="File containing target tokens")
    p.add_argument("--outgroup-token", action="append", default=[], help="Outgroup token excluded from event calling")
    p.add_argument("--outgroup-tokens-file", help="File containing outgroup tokens")
    p.add_argument("--exclude-token", action="append", default=[], help="Token excluded from event calling")
    p.add_argument("--exclude-tokens-file", help="File containing excluded tokens")
    p.add_argument("--target-state-mode", choices=["uniform", "allow-diverse"], default="uniform", help="Deprecated compatibility option; v0.1.5 reports identical/divergent target states automatically")
    p.add_argument("--min-target-non-gap", default="all", help="Deprecated alias for --min-target-informative")
    p.add_argument("--min-target-informative", help="Minimum informative target states; valid codons and full codon gaps are informative")
    p.add_argument("--max-target-gap-fraction-for-substitution", type=float, default=0.0)
    p.add_argument("--min-background-non-gap", type=int, default=5, help="Deprecated alias for --min-background-informative")
    p.add_argument("--min-background-informative", help="Minimum informative background states; valid codons and full codon gaps are informative")
    p.add_argument("--min-background-gap-fraction-for-target-non-gap-event", type=float, default=0.8)
    p.add_argument("--min-background-non-gap-fraction-for-target-gap-event", type=float, default=0.8)
    p.add_argument("--min-target-non-gap-fraction", type=float, default=1.0)
    p.add_argument("--min-target-gap-fraction", type=float, default=1.0)
    p.add_argument("--bed-mode", choices=["auto", "all-coordinateable", "substitution-only", "none"], default="auto")
    p.add_argument("--codon-table", choices=["universal", "vmitochondria"], default="universal")
    p.add_argument("--symbols", nargs="*", help="Only process selected symbols/stems; comma-separated values are accepted")
    p.add_argument("--symbols-file", help="File containing symbols to process")
    p.add_argument("--limit", type=int)
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_variants)

    p = sub.add_parser("verify", help="Verify code plus run, matrix, alignment, and variant outputs")
    p.add_argument(
        "--manifest",
        help="Species manifest TSV used to check FASTA headers; defaults to ./config/species_manifest.tsv if present",
    )
    p.add_argument("--output-root", help="Run output root containing fastas/, reports/, and optional matrices; defaults to current directory")
    p.add_argument("--full", action="store_true", help="Check every FASTA instead of a 200-file sample")
    p.add_argument(
        "--matrix-rows",
        choices=["none", "sample", "full"],
        default="sample",
        help="Matrix verification depth: none skips matrices, sample checks selected genes, full reads every matrix",
    )
    p.add_argument("--alignment-dir", help="Optional align output directory to verify; relative paths are resolved under --output-root")
    p.add_argument("--variant-dir", help="Optional variants output directory to verify; relative paths are resolved under --output-root")
    p.set_defaults(func=cmd_verify)

    p = sub.add_parser("summary", help="Print generated summary JSON files")
    p.add_argument("--output-root", help="Output root containing reports/; defaults to current directory")
    p.set_defaults(func=cmd_summary)

    p = sub.add_parser("test", help="Run packaged examples for reviewer-friendly smoke tests")
    p.add_argument("--example", choices=["mini"], default="mini")
    p.add_argument("--output-root", help="Where to write the mini run output; defaults to a temporary directory")
    p.add_argument("--force", action="store_true", help="Remove/rebuild output-root when supplied")
    p.set_defaults(func=cmd_test)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
