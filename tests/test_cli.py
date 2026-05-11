import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "refseq2cds.py"


def run_cli(*args: str, cwd: Path = ROOT) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(CLI), *args],
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    )


def test_help_lists_core_commands() -> None:
    proc = run_cli("--help")
    assert "manifest-from-gcf" in proc.stdout
    assert "run" in proc.stdout
    assert "test" in proc.stdout


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
