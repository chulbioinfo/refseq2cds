import csv
import gzip
import json
import importlib.util
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "refseq2cds.py"

_SPEC = importlib.util.spec_from_file_location("refseq2cds_cli", CLI)
assert _SPEC is not None and _SPEC.loader is not None
refseq2cds_cli = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(refseq2cds_cli)

_ALIGN_SPEC = importlib.util.spec_from_file_location(
    "refseq2cds_align",
    ROOT / "workflow" / "scripts" / "align_mafft_pal2nal.py",
)
assert _ALIGN_SPEC is not None and _ALIGN_SPEC.loader is not None
refseq2cds_align = importlib.util.module_from_spec(_ALIGN_SPEC)
sys.modules[_ALIGN_SPEC.name] = refseq2cds_align
_ALIGN_SPEC.loader.exec_module(refseq2cds_align)


def run_cli(*args: str, cwd: Path = ROOT) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(CLI), *args],
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    )


def fasta_headers(path: Path) -> list[str]:
    headers: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                headers.append(line[1:].strip().split()[0])
    return headers


def test_help_lists_core_commands() -> None:
    proc = run_cli("--help")
    assert "manifest-from-gcf" in proc.stdout
    assert "run" in proc.stdout
    assert "variants" in proc.stdout
    assert "test" in proc.stdout


def test_version_reports_current_version() -> None:
    proc = run_cli("--version")
    assert proc.stdout.strip() == "refseq2cds 0.1.5"


def test_gcf_accession_normalization_accepts_release_suffixes() -> None:
    assert refseq2cds_cli.normalize_gcf_accession("GCF_037993035.2-RS_2025_03") == "GCF_037993035.2"
    assert refseq2cds_cli.normalize_gcf_accession("GCF_964374335.1/") == "GCF_964374335.1"


def test_align_fasta_discovery_ignores_macos_appledouble_files(tmp_path: Path) -> None:
    (tmp_path / "GENE.fasta").write_text(">human\nATG\n")
    (tmp_path / "._GENE.fasta").write_text("not a fasta\n")

    selected = refseq2cds_align.collect_fastas(tmp_path, None, None)

    assert [path.name for path in selected] == ["GENE.fasta"]


def test_mini_example_generates_fastas_and_matrices(tmp_path: Path) -> None:
    out = tmp_path / "mini_run"
    proc = run_cli("test", "--example", "mini", "--output-root", str(out), "--force")
    result = json.loads(proc.stdout[proc.stdout.index("{") :])
    assert result["status"] == "pass"
    assert result["fastas"] == 2
    assert result["human_matrices"] == 2

    fasta_names = sorted(path.name for path in (out / "fastas").glob("*.fasta"))
    assert fasta_names == ["NEG1.fasta", "PASS1.fasta"]

    summary = json.loads((out / "reports" / "summary.json").read_text())
    assert summary["cds_qc_fail_reasons"] == {"internal_stop": 1}
    assert summary["rejected_reasons"]["component_size_not_expected_species_count"] == 2


def test_verify_can_skip_matrices_for_fasta_only_run(tmp_path: Path) -> None:
    out = tmp_path / "mini_no_matrices"
    run_cli(
        "run",
        "--manifest",
        str(ROOT / "examples/mini/manifest.tsv"),
        "--input-root",
        str(ROOT / "examples/mini/inputs"),
        "--output-root",
        str(out),
        "--offline",
        "--steps",
        "all",
        "--force",
    )

    proc = run_cli(
        "verify",
        "--output-root",
        str(out),
        "--manifest",
        str(ROOT / "examples/mini/manifest.tsv"),
        "--matrix-rows",
        "none",
    )
    result = json.loads(proc.stdout[proc.stdout.index("{") :])
    assert result["status"] == "pass"
    assert result["matrices"]["status"] == "skipped"


def test_mini_negative_strand_matrix_coordinates(tmp_path: Path) -> None:
    out = tmp_path / "mini_run"
    run_cli("test", "--example", "mini", "--output-root", str(out), "--force")

    matrix = out / "human_cds_matrices" / "NEG1.human_cds_genomic_matrix.tsv.gz"
    with gzip.open(matrix, "rt", newline="") as fh:
        first = next(csv.DictReader(fh, delimiter="\t"))

    assert first["strand"] == "-"
    assert first["cds_nt_pos_1based"] == "1"
    assert first["genomic_pos_1based"] == "562"
    assert first["bed_start_0based"] == "561"
    assert first["bed_end_0based"] == "562"


def test_reference_gene_mode_keeps_present_species_when_target_missing(tmp_path: Path) -> None:
    out = tmp_path / "reference_missing"
    run_cli(
        "run",
        "--manifest",
        str(ROOT / "examples/mini/manifest.tsv"),
        "--input-root",
        str(ROOT / "examples/mini/inputs"),
        "--output-root",
        str(out),
        "--offline",
        "--steps",
        "all",
        "--orthology-mode",
        "reference_gene_1to1_present_species",
        "--reference-symbol",
        "MISSING1",
        "--force",
    )

    fasta = out / "fastas" / "MISSING1.reference_1to1.cds.fasta"
    assert fasta_headers(fasta) == ["human", "chimpanzee"]

    query_dir = out / "results" / "reference_gene_1to1" / "MISSING1"
    assert (query_dir / "reference_gene.tsv").exists()
    assert (query_dir / "MISSING1.reference_1to1.meta.tsv").exists()
    rejected = (query_dir / "ortholog_rejected_graph.tsv").read_text()
    assert "gorilla" in rejected
    assert "no_ortholog_edge" in rejected


def test_reference_gene_mode_rejects_one_to_many_target_species(tmp_path: Path) -> None:
    out = tmp_path / "reference_paralog"
    run_cli(
        "run",
        "--manifest",
        str(ROOT / "examples/mini/manifest.tsv"),
        "--input-root",
        str(ROOT / "examples/mini/inputs"),
        "--output-root",
        str(out),
        "--offline",
        "--steps",
        "all",
        "--orthology-mode",
        "reference_gene_1to1_present_species",
        "--reference-symbol",
        "PARALOG1",
        "--force",
    )

    fasta = out / "fastas" / "PARALOG1.reference_1to1.cds.fasta"
    assert fasta_headers(fasta) == ["human", "chimpanzee"]

    query_dir = out / "results" / "reference_gene_1to1" / "PARALOG1"
    rejected = (query_dir / "ortholog_rejected_graph.tsv").read_text()
    assert "gorilla" in rejected
    assert "one_reference_to_many_targets" in rejected


def test_reference_gene_mode_passes_full_1to1_query(tmp_path: Path) -> None:
    out = tmp_path / "reference_pass"
    run_cli(
        "run",
        "--manifest",
        str(ROOT / "examples/mini/manifest.tsv"),
        "--input-root",
        str(ROOT / "examples/mini/inputs"),
        "--output-root",
        str(out),
        "--offline",
        "--steps",
        "all",
        "--orthology-mode",
        "reference_gene_1to1_present_species",
        "--reference-symbol",
        "PASS1",
        "--force",
    )

    fasta = out / "fastas" / "PASS1.reference_1to1.cds.fasta"
    assert fasta_headers(fasta) == ["human", "chimpanzee", "gorilla"]

    summary = json.loads((out / "results" / "reference_gene_1to1" / "PASS1" / "summary.json").read_text())
    assert summary["orthology_mode"] == "reference_gene_1to1_present_species"
    assert summary["final_fasta_sequence_count"] == 3
