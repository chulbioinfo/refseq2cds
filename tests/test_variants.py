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


def write_alignment(path: Path, records: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for token, seq in records.items():
            fh.write(f">{token}\n{seq}\n")


def write_human_map(path: Path, human_seq: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cds_codon = 0
    with gzip.open(path, "wt", newline="") as fh:
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
        for aln_idx, i in enumerate(range(0, len(human_seq), 3), start=1):
            codon = human_seq[i : i + 3]
            if codon == "---":
                writer.writerow(
                    {
                        "family_id": "GENE",
                        "token": "human",
                        "aln_codon_index_1based": aln_idx,
                        "cds_codon_index_1based": "",
                        "cds_nt_start_1based": "",
                        "cds_nt_end_1based": "",
                        "is_gap": "true",
                    }
                )
            else:
                cds_codon += 1
                writer.writerow(
                    {
                        "family_id": "GENE",
                        "token": "human",
                        "aln_codon_index_1based": aln_idx,
                        "cds_codon_index_1based": cds_codon,
                        "cds_nt_start_1based": (cds_codon - 1) * 3 + 1,
                        "cds_nt_end_1based": cds_codon * 3,
                        "is_gap": "false",
                    }
                )


def write_human_matrix(path: Path, nt_count: int = 12, strand: str = "+") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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
    ]
    with gzip.open(path, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for i in range(1, nt_count + 1):
            writer.writerow(
                {
                    "family_id": "GENE",
                    "human_symbol": "GENE",
                    "human_GeneID": "1",
                    "transcript_accession": "NM_TEST",
                    "protein_accession": "NP_TEST",
                    "cds_nt_pos_1based": i,
                    "codon_index_1based": ((i - 1) // 3) + 1,
                    "codon_offset_1based": ((i - 1) % 3) + 1,
                    "cds_base_transcript_orientation": "A",
                    "refseq_accession": "NC_TEST",
                    "ucsc_chrom": "chrTest",
                    "genomic_pos_1based": 1000 + i,
                    "bed_start_0based": 999 + i,
                    "bed_end_0based": 1000 + i,
                    "strand": strand,
                }
            )


def prepare_variant_fixture(tmp_path: Path, records: dict[str, str]) -> tuple[Path, Path, Path, Path]:
    aln = tmp_path / "alignments" / "GENE.codon.fasta"
    maps = tmp_path / "maps" / "GENE.human.codon_map.tsv.gz"
    matrices = tmp_path / "matrices" / "GENE.human_cds_genomic_matrix.tsv.gz"
    out = tmp_path / "variants"
    write_alignment(aln, records)
    write_human_map(maps, records["human"])
    write_human_matrix(matrices)
    return aln.parent, maps.parent, matrices.parent, out


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_variants_default_uniform_calls_human_specific_substitution(tmp_path: Path) -> None:
    aln_dir, map_dir, matrix_dir, out = prepare_variant_fixture(
        tmp_path,
        {
            "human": "AAA",
            "chimpanzee": "CAA",
            "gorilla": "GAA",
            "out": "AAA",
        },
    )
    (aln_dir / "._GENE.codon.fasta").write_text("not a fasta\n")
    proc = run_cli(
        "variants",
        "--alignment-dir",
        str(aln_dir),
        "--codon-map-dir",
        str(map_dir),
        "--matrix-dir",
        str(matrix_dir),
        "--output-dir",
        str(out),
        "--coordinate-reference-token",
        "human",
        "--target-token",
        "human",
        "--outgroup-token",
        "out",
        "--min-background-non-gap",
        "2",
        "--force",
    )
    summary = json.loads(proc.stdout[proc.stdout.index("{") :])
    assert summary["target_state_mode"] == "uniform"
    assert summary["variant_count"] == 1
    assert summary["identical_sequence_count"] == 1
    assert summary["bed_rows"] == 1
    assert summary["merged_bed_path"] == "merged.bed"
    assert summary["merged_bed_rows"] == 1
    rows = read_gzip_tsv(out / "matrices" / "GENE.variant_matrix.tsv.gz")
    assert rows[0]["coordinateable"] == "true"
    assert rows[0]["event_type"] == "aa_variant"
    assert rows[0]["event_subtype"] == "identical_sequence"
    assert rows[0]["bed_event_class"] == "variant"
    assert (out / "merged.bed").read_text().startswith("chrTest\t1000\t1003\tGENE|target=human|variant")


def test_variants_default_calls_divergent_target_variant(tmp_path: Path) -> None:
    records = {
        "human": "AAA",
        "chimpanzee": "AGA",
        "gorilla": "CAA",
        "bonobo": "GAA",
    }
    aln_dir, map_dir, matrix_dir, out = prepare_variant_fixture(tmp_path, records)
    proc = run_cli(
        "variants",
        "--alignment-dir",
        str(aln_dir),
        "--codon-map-dir",
        str(map_dir),
        "--matrix-dir",
        str(matrix_dir),
        "--output-dir",
        str(out),
        "--coordinate-reference-token",
        "human",
        "--target-token",
        "human",
        "--target-token",
        "chimpanzee",
        "--min-background-non-gap",
        "2",
        "--force",
    )
    summary = json.loads(proc.stdout[proc.stdout.index("{") :])
    assert summary["variant_count"] == 1
    assert summary["divergent_sequence_count"] == 1
    rows = read_gzip_tsv(out / "matrices" / "GENE.variant_matrix.tsv.gz")
    assert rows[0]["event_subtype"] == "divergent_sequence"


def test_variants_rejects_non_exclusive_gap_rich_background(tmp_path: Path) -> None:
    aln_dir, map_dir, matrix_dir, out = prepare_variant_fixture(
        tmp_path,
        {
            "human": "AAA",
            "chimpanzee": "AAA",
            "gorilla": "---",
            "bonobo": "---",
        },
    )
    proc = run_cli(
        "variants",
        "--alignment-dir",
        str(aln_dir),
        "--codon-map-dir",
        str(map_dir),
        "--matrix-dir",
        str(matrix_dir),
        "--output-dir",
        str(out),
        "--coordinate-reference-token",
        "human",
        "--target-token",
        "human",
        "--min-background-informative",
        "1",
        "--force",
    )
    summary = json.loads(proc.stdout[proc.stdout.index("{") :])
    assert summary["variant_count"] == 0
    assert summary["bed_rows"] == 0
    assert summary["merged_bed_rows"] == 0
    assert (out / "merged.bed").read_text() == ""


def test_variants_coordinate_reference_gap_keeps_matrix_only_event(tmp_path: Path) -> None:
    aln_dir, map_dir, matrix_dir, out = prepare_variant_fixture(
        tmp_path,
        {
            "human": "---",
            "chimpanzee": "AAA",
            "gorilla": "CAA",
            "bonobo": "GAA",
        },
    )
    proc = run_cli(
        "variants",
        "--alignment-dir",
        str(aln_dir),
        "--codon-map-dir",
        str(map_dir),
        "--matrix-dir",
        str(matrix_dir),
        "--output-dir",
        str(out),
        "--coordinate-reference-token",
        "human",
        "--target-token",
        "chimpanzee",
        "--min-background-non-gap",
        "2",
        "--force",
    )
    summary = json.loads(proc.stdout[proc.stdout.index("{") :])
    assert summary["event_count"] >= 1
    assert summary["variant_count"] == 1
    assert summary["bed_rows"] == 0
    rows = read_gzip_tsv(out / "matrices" / "GENE.variant_matrix.tsv.gz")
    assert any(row["coordinate_status"] == "not_mapped_coordinate_reference_gap" for row in rows)


def test_variants_target_gap_maps_as_variant_when_reference_has_base(tmp_path: Path) -> None:
    aln_dir, map_dir, matrix_dir, out = prepare_variant_fixture(
        tmp_path,
        {
            "human": "AAA",
            "chimpanzee": "---",
            "gorilla": "AAA",
            "bonobo": "AAA",
        },
    )
    proc = run_cli(
        "variants",
        "--alignment-dir",
        str(aln_dir),
        "--codon-map-dir",
        str(map_dir),
        "--matrix-dir",
        str(matrix_dir),
        "--output-dir",
        str(out),
        "--coordinate-reference-token",
        "human",
        "--target-token",
        "chimpanzee",
        "--min-background-non-gap",
        "2",
        "--force",
    )
    summary = json.loads(proc.stdout[proc.stdout.index("{") :])
    assert summary["variant_count"] == 1
    assert summary["bed_rows"] == 1
    assert summary["merged_bed_rows"] == 1
    bed = out / "bed" / "GENE.chimpanzee.variant.bed"
    assert bed.read_text().startswith("chrTest\t1000\t1003\tGENE|target=chimpanzee|variant")
    assert (out / "merged.bed").read_text() == bed.read_text()
