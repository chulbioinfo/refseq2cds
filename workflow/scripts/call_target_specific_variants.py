#!/usr/bin/env python3
"""
Target-vs-background comparative CDS variant scanner.

This driver consumes codon-aware alignments, alignment-to-CDS codon maps, and a
coordinate-reference CDS-to-genome matrix. It calls target-group-specific amino
acid variants by requiring mutually exclusive target and background amino acid
state sets; the coordinate reference is used only for coordinate mapping and
BED output.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
import time
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from Bio import SeqIO
from Bio.Data.CodonTable import unambiguous_dna_by_id


ROOT = Path(__file__).resolve().parents[2]
STANDARD_TABLE = unambiguous_dna_by_id[1]
VERTEBRATE_MITO_TABLE = unambiguous_dna_by_id[2]


def log(message: str) -> None:
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def rel(path: Path, base: Path) -> str:
    try:
        return str(path.relative_to(base))
    except ValueError:
        return str(path)


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
    for record in SeqIO.parse(str(path), "fasta"):
        ident = record.id.split()[0]
        if ident in records:
            raise ValueError(f"Duplicate FASTA ID {ident} in {path}")
        records[ident] = str(record.seq).upper().replace("U", "T")
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    lengths = {len(seq) for seq in records.values()}
    if len(lengths) != 1:
        raise ValueError(f"Sequences in {path} do not have equal alignment lengths: {sorted(lengths)}")
    length = next(iter(lengths))
    if length % 3 != 0:
        raise ValueError(f"Codon alignment length is not a multiple of 3 in {path}: {length}")
    return records


def read_token_file(path: Optional[str]) -> List[str]:
    if not path:
        return []
    tokens: List[str] = []
    with Path(path).open() as fh:
        for line in fh:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            tokens.extend(part for part in text.replace(",", " ").split() if part)
    return tokens


def unique_preserve_order(values: Iterable[str]) -> List[str]:
    seen: Set[str] = set()
    out: List[str] = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return out


def split_values(values: Optional[Sequence[str]]) -> List[str]:
    if not values:
        return []
    out: List[str] = []
    for value in values:
        out.extend(part for part in value.replace(",", " ").split() if part)
    return out


def symbol_from_alignment(path: Path) -> str:
    name = path.name
    for suffix in [
        ".reference_1to1.cds.codon.fasta",
        ".reference_1to1.cds.fasta",
        ".codon.fasta",
        ".cds.fasta",
        ".fasta",
        ".fa",
    ]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def sanitize_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value).strip("_") or "target"


def collect_alignment_files(alignment_dir: Path, symbols: Sequence[str], limit: Optional[int]) -> List[Path]:
    if (alignment_dir / "codon").is_dir():
        alignment_dir = alignment_dir / "codon"
    candidates = sorted(
        p
        for p in alignment_dir.iterdir()
        if (
            p.is_file()
            and not p.name.startswith("._")
            and (p.name.endswith(".codon.fasta") or p.name.endswith(".fasta") or p.name.endswith(".fa"))
        )
    )
    wanted = set(symbols)
    if wanted:
        candidates = [p for p in candidates if symbol_from_alignment(p) in wanted or p.stem in wanted]
        observed = {symbol_from_alignment(p) for p in candidates} | {p.stem for p in candidates}
        missing = sorted(wanted - observed)
        if missing:
            raise FileNotFoundError(f"Requested symbols not found in {alignment_dir}: {missing}")
    if limit is not None:
        candidates = candidates[:limit]
    if not candidates:
        raise FileNotFoundError(f"No codon alignment files found in {alignment_dir}")
    return candidates


def codon_table(name: str):
    if name == "universal":
        return STANDARD_TABLE
    if name == "vmitochondria":
        return VERTEBRATE_MITO_TABLE
    raise ValueError(f"Unsupported codon table: {name}")


@dataclass(frozen=True)
class CodonState:
    token: str
    codon: str
    status: str
    aa: str = ""

    @property
    def is_valid(self) -> bool:
        return self.status == "valid_codon"

    @property
    def is_gap(self) -> bool:
        return self.status == "full_gap"

    @property
    def is_informative_aa(self) -> bool:
        return self.status in {"valid_codon", "full_gap"}


def classify_codon(token: str, codon: str, table) -> CodonState:
    codon = codon.upper().replace("U", "T")
    if len(codon) != 3:
        return CodonState(token, codon, "invalid_length")
    if codon == "---":
        return CodonState(token, codon, "full_gap", "-")
    if "-" in codon:
        return CodonState(token, codon, "partial_gap")
    if any(base not in {"A", "C", "G", "T"} for base in codon):
        return CodonState(token, codon, "ambiguous_nt")
    if codon in set(table.stop_codons):
        return CodonState(token, codon, "stop_codon", "*")
    aa = table.forward_table.get(codon, "X")
    if aa == "X":
        return CodonState(token, codon, "ambiguous_nt", aa)
    return CodonState(token, codon, "valid_codon", aa)


@dataclass
class GroupStats:
    tokens: List[str]
    states: List[CodonState]

    @property
    def count(self) -> int:
        return len(self.states)

    @property
    def valid(self) -> List[CodonState]:
        return [s for s in self.states if s.is_valid]

    @property
    def gaps(self) -> List[CodonState]:
        return [s for s in self.states if s.is_gap]

    @property
    def ambiguous(self) -> List[CodonState]:
        return [s for s in self.states if s.status in {"partial_gap", "ambiguous_nt", "stop_codon", "invalid_length"}]

    @property
    def aa_set(self) -> Set[str]:
        return {s.aa for s in self.informative}

    @property
    def codon_set(self) -> Set[str]:
        return {s.codon for s in self.informative}

    @property
    def informative(self) -> List[CodonState]:
        return [s for s in self.states if s.is_informative_aa]

    def nt_set(self, offset_1based: int) -> Set[str]:
        return {s.codon[offset_1based - 1] for s in self.informative}

    @property
    def gap_fraction(self) -> float:
        return len(self.gaps) / self.count if self.count else 0.0

    @property
    def non_gap_fraction(self) -> float:
        return len(self.valid) / self.count if self.count else 0.0


@dataclass
class TokenGroups:
    all_tokens: List[str]
    target_tokens: List[str]
    background_tokens: List[str]
    outgroup_tokens: List[str]
    exclude_tokens: List[str]
    coordinate_reference_token: str
    coordinate_reference_role: str


def resolve_token_groups(
    all_tokens: Sequence[str],
    target_tokens: Sequence[str],
    outgroup_tokens: Sequence[str],
    exclude_tokens: Sequence[str],
    coordinate_reference_token: str,
) -> TokenGroups:
    all_set = set(all_tokens)
    targets = unique_preserve_order(target_tokens)
    outgroups = unique_preserve_order(outgroup_tokens)
    excludes = unique_preserve_order(exclude_tokens)
    if not targets:
        raise ValueError("At least one --target-token or --target-tokens-file entry is required")
    missing_targets = sorted(set(targets) - all_set)
    if missing_targets:
        raise ValueError(f"Target tokens are missing from alignment: {missing_targets}")
    overlaps = {
        "target/outgroup": set(targets) & set(outgroups),
        "target/exclude": set(targets) & set(excludes),
        "outgroup/exclude": set(outgroups) & set(excludes),
    }
    bad = {k: sorted(v) for k, v in overlaps.items() if v}
    if bad:
        raise ValueError(f"Token groups must be disjoint: {bad}")
    background = [t for t in all_tokens if t not in set(targets) | set(outgroups) | set(excludes)]
    if not background:
        raise ValueError("Background token set is empty after removing target/outgroup/exclude tokens")
    if coordinate_reference_token in targets:
        role = "target"
    elif coordinate_reference_token in background:
        role = "background"
    elif coordinate_reference_token in outgroups:
        role = "outgroup"
    elif coordinate_reference_token in excludes:
        role = "excluded"
    else:
        role = "absent"
    return TokenGroups(
        all_tokens=list(all_tokens),
        target_tokens=targets,
        background_tokens=background,
        outgroup_tokens=outgroups,
        exclude_tokens=excludes,
        coordinate_reference_token=coordinate_reference_token,
        coordinate_reference_role=role,
    )


def threshold_count(value: str, total: int, default_all: bool = False) -> int:
    text = str(value).strip().lower()
    if text == "all":
        return total
    if not text:
        return total if default_all else 0
    try:
        if "." in text:
            frac = float(text)
            if 0 <= frac <= 1:
                return int(total * frac + 0.999999)
        return int(text)
    except ValueError as exc:
        raise ValueError(f"Invalid threshold count: {value}") from exc


@dataclass
class Event:
    event_id: str
    family_id: str
    symbol: str
    source_alignment_file: str
    event_type: str
    event_subtype: str
    target_state_mode: str
    bed_mode: str
    aln_codon_start_1based: int
    aln_codon_end_1based: int
    target_stats: GroupStats
    background_stats: GroupStats
    outgroup_stats: GroupStats
    coordinate_ref_state: Optional[CodonState]
    target_state_set: Set[str]
    background_state_set: Set[str]
    target_codon_state_set: Set[str]
    background_codon_state_set: Set[str]
    event_call_basis: str = "target_vs_background_exclusive"
    is_nonsynonymous: bool = False
    codon_change_class: str = ""
    aa_change_label: str = ""
    polarity_inferred: bool = False
    ancestral_state_used: bool = False


@dataclass
class NtChange:
    event: Event
    codon_offset_1based: int
    target_nt_state_set: Set[str]
    background_nt_state_set: Set[str]
    nt_change_role: str


@dataclass
class MappingResources:
    codon_map: Dict[int, dict] = field(default_factory=dict)
    matrix_by_nt: Dict[int, dict] = field(default_factory=dict)
    codon_map_path: str = ""
    matrix_path: str = ""


def load_codon_map(path: Path) -> Dict[int, dict]:
    out: Dict[int, dict] = {}
    with gzip.open(path, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[int(row["aln_codon_index_1based"])] = row
    return out


def load_matrix(path: Path) -> Dict[int, dict]:
    out: Dict[int, dict] = {}
    with gzip.open(path, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[int(row["cds_nt_pos_1based"])] = row
    return out


def find_codon_map(codon_map_dir: Path, symbol: str, alignment_path: Path, token: str) -> Optional[Path]:
    if not codon_map_dir.exists():
        return None
    stems = [alignment_path.stem]
    if stems[0].endswith(".codon"):
        stems.append(stems[0][: -len(".codon")])
    stems.append(symbol)
    seen = set()
    for stem in stems:
        if stem in seen:
            continue
        seen.add(stem)
        candidate = codon_map_dir / f"{stem}.{token}.codon_map.tsv.gz"
        if candidate.exists():
            return candidate
    matches = sorted(codon_map_dir.glob(f"*{symbol}*.{token}.codon_map.tsv.gz"))
    return matches[0] if matches else None


def find_matrix(matrix_dir: Path, symbol: str, token: str) -> Optional[Path]:
    if not matrix_dir.exists():
        return None
    candidates = [
        matrix_dir / f"{symbol}.{token}_cds_genomic_matrix.tsv.gz",
        matrix_dir / f"{symbol}.human_cds_genomic_matrix.tsv.gz",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    matches = sorted(matrix_dir.glob(f"{symbol}.*cds_genomic_matrix.tsv.gz"))
    return matches[0] if matches else None


def load_mapping_resources(
    codon_map_dir: Path,
    matrix_dir: Path,
    symbol: str,
    alignment_path: Path,
    token: str,
    output_root: Path,
) -> MappingResources:
    resources = MappingResources()
    map_path = find_codon_map(codon_map_dir, symbol, alignment_path, token)
    matrix_path = find_matrix(matrix_dir, symbol, token)
    if map_path:
        resources.codon_map = load_codon_map(map_path)
        resources.codon_map_path = rel(map_path, output_root)
    if matrix_path:
        resources.matrix_by_nt = load_matrix(matrix_path)
        resources.matrix_path = rel(matrix_path, output_root)
    return resources


def state_csv(values: Iterable[str]) -> str:
    return ",".join(sorted(str(v) for v in values if v != ""))


def tokens_csv(values: Sequence[str]) -> str:
    return ",".join(values)


def group_stats_rows(stats: GroupStats, prefix: str) -> Dict[str, str]:
    return {
        f"{prefix}_count": str(stats.count),
        f"{prefix}_valid_count": str(len(stats.valid)),
        f"{prefix}_informative_count": str(len(stats.informative)),
        f"{prefix}_gap_count": str(len(stats.gaps)),
        f"{prefix}_ambiguous_count": str(len(stats.ambiguous)),
        f"{prefix}_state_set": state_csv(stats.aa_set),
        f"{prefix}_codon_state_set": state_csv(stats.codon_set),
        f"{prefix}_uniform": str(len(stats.aa_set) == 1).lower(),
    }


def event_base_row(event: Event, groups: TokenGroups) -> Dict[str, str]:
    coord = event.coordinate_ref_state
    row = {
        "family_id": event.family_id,
        "reference_symbol": event.symbol,
        "source_alignment_file": event.source_alignment_file,
        "event_id": event.event_id,
        "event_type": event.event_type,
        "event_subtype": event.event_subtype,
        "target_state_mode": event.target_state_mode,
        "bed_mode": event.bed_mode,
        "all_tokens": tokens_csv(groups.all_tokens),
        "target_tokens": tokens_csv(groups.target_tokens),
        "background_tokens": tokens_csv(groups.background_tokens),
        "outgroup_tokens": tokens_csv(groups.outgroup_tokens),
        "exclude_tokens": tokens_csv(groups.exclude_tokens),
        "coordinate_reference_token": groups.coordinate_reference_token,
        "coordinate_reference_role": groups.coordinate_reference_role,
        "aln_codon_start_1based": str(event.aln_codon_start_1based),
        "aln_codon_end_1based": str(event.aln_codon_end_1based),
        "aln_nt_start_1based": str((event.aln_codon_start_1based - 1) * 3 + 1),
        "aln_nt_end_1based": str(event.aln_codon_end_1based * 3),
        "target_state_set": state_csv(event.target_state_set),
        "target_codon_state_set": state_csv(event.target_codon_state_set),
        "background_state_set": state_csv(event.background_state_set),
        "background_codon_state_set": state_csv(event.background_codon_state_set),
        "background_contains_target_state": str(bool(event.target_state_set & event.background_state_set)).lower(),
        "outgroup_state_set": state_csv(event.outgroup_stats.aa_set),
        "event_call_basis": event.event_call_basis,
        "is_nonsynonymous": str(event.is_nonsynonymous).lower(),
        "codon_change_class": event.codon_change_class,
        "nt_change_role": "",
        "aa_change_label": event.aa_change_label,
        "polarity_inferred": str(event.polarity_inferred).lower(),
        "ancestral_state_used": str(event.ancestral_state_used).lower(),
        "coordinate_ref_codon": coord.codon if coord else "",
        "coordinate_ref_aa": coord.aa if coord else "",
        "coordinate_ref_codon_status": coord.status if coord else "missing_token",
        "coordinate_ref_nt_cds_strand": "",
        "coordinate_ref_nt_genome_plus": "",
        "coordinate_ref_has_mappable_base": "false",
        "coordinate_ref_cds_codon_start_1based": "",
        "coordinate_ref_cds_codon_end_1based": "",
        "coordinate_ref_cds_nt_start_1based": "",
        "coordinate_ref_cds_nt_end_1based": "",
        "codon_offset_1based": "",
        "target_nt_state_set": "",
        "background_nt_state_set": "",
        "coordinateable": "false",
        "coordinate_status": "",
        "ucsc_chrom": "",
        "bed_start_0based": "",
        "bed_end_0based": "",
        "strand": "",
        "bed_event_class": "",
        "event_block_count": "",
        "event_block_index": "",
        "genomic_block_id": "",
    }
    row.update(group_stats_rows(event.target_stats, "target"))
    row.update(group_stats_rows(event.background_stats, "background"))
    return row


AA_EVENT_COLUMNS = [
    "family_id",
    "reference_symbol",
    "source_alignment_file",
    "event_id",
    "event_type",
    "event_subtype",
    "target_state_mode",
    "bed_mode",
    "all_tokens",
    "target_tokens",
    "background_tokens",
    "outgroup_tokens",
    "exclude_tokens",
    "coordinate_reference_token",
    "coordinate_reference_role",
    "aln_codon_start_1based",
    "aln_codon_end_1based",
    "aln_nt_start_1based",
    "aln_nt_end_1based",
    "target_count",
    "target_valid_count",
    "target_informative_count",
    "target_gap_count",
    "target_ambiguous_count",
    "target_state_set",
    "target_codon_state_set",
    "target_uniform",
    "background_count",
    "background_valid_count",
    "background_informative_count",
    "background_gap_count",
    "background_ambiguous_count",
    "background_state_set",
    "background_codon_state_set",
    "background_uniform",
    "background_contains_target_state",
    "outgroup_state_set",
    "coordinate_ref_codon",
    "coordinate_ref_aa",
    "coordinate_ref_codon_status",
    "event_call_basis",
    "is_nonsynonymous",
    "codon_change_class",
    "aa_change_label",
    "polarity_inferred",
    "ancestral_state_used",
]

MATRIX_COLUMNS = AA_EVENT_COLUMNS + [
    "nt_change_role",
    "coordinate_ref_nt_cds_strand",
    "coordinate_ref_nt_genome_plus",
    "coordinate_ref_has_mappable_base",
    "coordinate_ref_cds_codon_start_1based",
    "coordinate_ref_cds_codon_end_1based",
    "coordinate_ref_cds_nt_start_1based",
    "coordinate_ref_cds_nt_end_1based",
    "codon_offset_1based",
    "target_nt_state_set",
    "background_nt_state_set",
    "coordinateable",
    "coordinate_status",
    "ucsc_chrom",
    "bed_start_0based",
    "bed_end_0based",
    "strand",
    "bed_event_class",
    "event_block_count",
    "event_block_index",
    "genomic_block_id",
]

NT_COLUMNS = [
    "family_id",
    "reference_symbol",
    "event_id",
    "event_type",
    "event_subtype",
    "target_state_mode",
    "aln_codon_index_1based",
    "codon_offset_1based",
    "target_nt_state_set",
    "background_nt_state_set",
    "codon_change_class",
    "nt_change_role",
    "is_nonsynonymous",
]


def write_gzip_tsv(path: Path, columns: Sequence[str], rows: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(columns), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def group_stats_for(states: Dict[str, CodonState], tokens: Sequence[str]) -> GroupStats:
    return GroupStats(tokens=list(tokens), states=[states[t] for t in tokens if t in states])


def make_event_id(n: int) -> str:
    return f"event_{n:06d}"


def call_aa_variant(
    event_id: str,
    family_id: str,
    symbol: str,
    source_alignment_file: str,
    aln_codon_index: int,
    target_stats: GroupStats,
    background_stats: GroupStats,
    outgroup_stats: GroupStats,
    coord_state: Optional[CodonState],
    target_state_mode: str,
    bed_mode: str,
    min_target_informative: int,
    min_background_informative: int,
) -> Optional[Event]:
    if len(target_stats.informative) < min_target_informative:
        return None
    if len(background_stats.informative) < min_background_informative:
        return None
    if not target_stats.aa_set:
        return None
    if not background_stats.aa_set:
        return None
    if target_stats.aa_set & background_stats.aa_set:
        return None
    subtype = "identical_sequence" if len(target_stats.aa_set) == 1 else "divergent_sequence"
    target_aa = state_csv(target_stats.aa_set)
    bg_aa = state_csv(background_stats.aa_set)
    codon_class = "aa_mutually_exclusive_variant"
    is_nonsynonymous = "-" not in target_stats.aa_set and "-" not in background_stats.aa_set
    return Event(
        event_id=event_id,
        family_id=family_id,
        symbol=symbol,
        source_alignment_file=source_alignment_file,
        event_type="aa_variant",
        event_subtype=subtype,
        target_state_mode=target_state_mode,
        bed_mode=bed_mode,
        aln_codon_start_1based=aln_codon_index,
        aln_codon_end_1based=aln_codon_index,
        target_stats=target_stats,
        background_stats=background_stats,
        outgroup_stats=outgroup_stats,
        coordinate_ref_state=coord_state,
        target_state_set=set(target_stats.aa_set),
        background_state_set=set(background_stats.aa_set),
        target_codon_state_set=set(target_stats.codon_set),
        background_codon_state_set=set(background_stats.codon_set),
        event_call_basis="target_vs_background_mutually_exclusive_aa_state",
        is_nonsynonymous=is_nonsynonymous,
        codon_change_class=codon_class,
        aa_change_label=f"target:{target_aa}|background:{bg_aa}",
    )


def complement_base(base: str) -> str:
    return base.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))


def matrix_nt_for_row(row: dict) -> Tuple[str, str]:
    cds_base = row.get("cds_base_transcript_orientation", "")
    strand = row.get("strand", "")
    if strand == "-":
        return cds_base, complement_base(cds_base)
    return cds_base, cds_base


def coordinate_for_nt(
    aln_codon_index: int,
    offset: int,
    event: Event,
    resources: MappingResources,
) -> Tuple[Optional[dict], str, str, str, str, str]:
    coord_state = event.coordinate_ref_state
    if coord_state is None:
        return None, "not_mapped_coordinate_reference_absent", "", "", "", ""
    if coord_state.status == "full_gap":
        return None, "not_mapped_coordinate_reference_gap", "", "", "", ""
    if coord_state.status == "partial_gap":
        return None, "not_mapped_coordinate_reference_partial_gap", "", "", "", ""
    if coord_state.status in {"ambiguous_nt", "stop_codon", "invalid_length"}:
        return None, f"not_mapped_coordinate_reference_{coord_state.status}", "", "", "", ""
    if not resources.codon_map:
        return None, "not_mapped_missing_codon_map", "", "", "", ""
    if not resources.matrix_by_nt:
        return None, "not_mapped_missing_reference_matrix", "", "", "", ""
    map_row = resources.codon_map.get(aln_codon_index)
    if not map_row:
        return None, "not_mapped_missing_codon_map", "", "", "", ""
    if map_row.get("is_gap", "").lower() == "true" or not map_row.get("cds_nt_start_1based"):
        return None, "not_mapped_coordinate_reference_gap", "", "", "", ""
    cds_nt = int(map_row["cds_nt_start_1based"]) + offset - 1
    matrix_row = resources.matrix_by_nt.get(cds_nt)
    if not matrix_row:
        return None, "not_mapped_no_genomic_coordinate", str(cds_nt), "", "", ""
    cds_base, genome_base = matrix_nt_for_row(matrix_row)
    return (
        matrix_row,
        "mapped_coordinate_reference_base",
        str(cds_nt),
        str(map_row.get("cds_codon_index_1based", "")),
        cds_base,
        genome_base,
    )


def contiguous_blocks(rows: List[dict]) -> List[List[dict]]:
    if not rows:
        return []
    sorted_rows = sorted(rows, key=lambda r: (r.get("ucsc_chrom", ""), int(r["bed_start_0based"])))
    blocks: List[List[dict]] = [[sorted_rows[0]]]
    for row in sorted_rows[1:]:
        prev = blocks[-1][-1]
        same = row.get("ucsc_chrom") == prev.get("ucsc_chrom") and row.get("strand") == prev.get("strand")
        adjacent = int(row["bed_start_0based"]) == int(prev["bed_end_0based"])
        if same and adjacent:
            blocks[-1].append(row)
        else:
            blocks.append([row])
    return blocks


def matrix_rows_for_event(
    event: Event,
    groups: TokenGroups,
    resources: MappingResources,
    output_root: Path,
) -> Tuple[List[dict], List[dict]]:
    base = event_base_row(event, groups)
    nt_rows: List[dict] = []
    matrix_rows: List[dict] = []
    mapped_nt_rows: List[dict] = []
    failure_status = ""
    for aln_idx in range(event.aln_codon_start_1based, event.aln_codon_end_1based + 1):
        for offset in [1, 2, 3]:
            matrix_row, status, _cds_nt, _cds_codon, _cds_base, _genome_base = coordinate_for_nt(
                aln_idx,
                offset,
                event,
                resources,
            )
            if matrix_row:
                mapped_nt_rows.append(matrix_row)
            elif not failure_status:
                failure_status = status
    if not mapped_nt_rows:
        row = dict(base)
        row["coordinate_status"] = failure_status or "not_mapped_no_genomic_coordinate"
        row["bed_event_class"] = "variant"
        matrix_rows.append(row)
        return matrix_rows, nt_rows

    blocks = contiguous_blocks(mapped_nt_rows)
    for i, block in enumerate(blocks, start=1):
        row = dict(base)
        row.update(
            {
                "coordinate_ref_has_mappable_base": "true",
                "coordinateable": "true",
                "coordinate_status": "mapped_coordinate_reference_base",
                "ucsc_chrom": block[0].get("ucsc_chrom", ""),
                "bed_start_0based": str(min(int(r["bed_start_0based"]) for r in block)),
                "bed_end_0based": str(max(int(r["bed_end_0based"]) for r in block)),
                "strand": block[0].get("strand", ""),
                "bed_event_class": "variant",
                "event_block_count": str(len(blocks)),
                "event_block_index": str(i),
                "genomic_block_id": f"{event.event_id}.block_{i}",
                "coordinate_ref_cds_nt_start_1based": block[0].get("cds_nt_pos_1based", ""),
                "coordinate_ref_cds_nt_end_1based": block[-1].get("cds_nt_pos_1based", ""),
                "coordinate_ref_cds_codon_start_1based": str(event.aln_codon_start_1based),
                "coordinate_ref_cds_codon_end_1based": str(event.aln_codon_end_1based),
            }
        )
        matrix_rows.append(row)
    return matrix_rows, nt_rows


def bed_rows_from_matrix(matrix_rows: Sequence[dict], symbol: str, target_label: str, bed_mode: str) -> Dict[str, List[List[str]]]:
    beds = {"variant": []}
    if bed_mode == "none":
        return beds
    for row in matrix_rows:
        if row.get("coordinateable") != "true":
            continue
        bed_class = row.get("bed_event_class", "")
        if bed_class not in beds:
            continue
        name = (
            f"{symbol}|target={target_label}|{bed_class}|{row.get('event_id')}|"
            f"block_{row.get('event_block_index')}_of_{row.get('event_block_count')}|"
            f"codon_{row.get('aln_codon_start_1based')}_{row.get('aln_codon_end_1based')}|"
            f"{row.get('event_subtype')}"
        )
        beds[bed_class].append(
            [
                row.get("ucsc_chrom", ""),
                row.get("bed_start_0based", ""),
                row.get("bed_end_0based", ""),
                name,
                "0",
                row.get("strand", ".") or ".",
            ]
        )
    return beds


def write_bed(path: Path, rows: Sequence[Sequence[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        for row in rows:
            if row[0] and row[1] and row[2] and int(row[1]) < int(row[2]):
                fh.write("\t".join(row) + "\n")


def bed_sort_key(row: Sequence[str]) -> Tuple[str, int, int, str]:
    return (row[0], int(row[1]), int(row[2]), row[3] if len(row) > 3 else "")


def write_merged_bed(output_root: Path, results: Sequence["GeneResult"]) -> Tuple[str, int]:
    rows: List[List[str]] = []
    for result in results:
        for path_text in result.bed_paths.split(";"):
            if not path_text:
                continue
            bed_path = output_root / path_text
            if not bed_path.exists():
                continue
            with bed_path.open(newline="") as fh:
                for line in fh:
                    text = line.rstrip("\n")
                    if not text:
                        continue
                    fields = text.split("\t")
                    if len(fields) >= 3:
                        rows.append(fields)
    merged_path = output_root / "merged.bed"
    write_bed(merged_path, sorted(rows, key=bed_sort_key))
    with merged_path.open() as fh:
        merged_rows = sum(1 for line in fh if line.strip())
    return rel(merged_path, output_root), merged_rows


@dataclass
class GeneResult:
    status: str
    symbol: str
    family_id: str
    alignment_path: str
    aa_events_path: str = ""
    codon_events_path: str = ""
    nt_changes_path: str = ""
    variant_matrix_path: str = ""
    bed_paths: str = ""
    event_count: int = 0
    variant_count: int = 0
    identical_sequence_count: int = 0
    divergent_sequence_count: int = 0
    substitution_count: int = 0
    indel_like_count: int = 0
    nt_change_count: int = 0
    coordinateable_rows: int = 0
    bed_rows: int = 0
    reason: str = ""


def process_gene(path: Path, args: argparse.Namespace, output_root: Path) -> GeneResult:
    symbol = symbol_from_alignment(path)
    records = read_fasta(path)
    all_tokens = list(records)
    groups = resolve_token_groups(
        all_tokens,
        args.target_tokens,
        args.outgroup_tokens,
        args.exclude_tokens,
        args.coordinate_reference_token,
    )
    table = codon_table(args.codon_table)
    resources = load_mapping_resources(
        Path(args.codon_map_dir),
        Path(args.matrix_dir),
        symbol,
        path,
        args.coordinate_reference_token,
        output_root,
    )
    source_alignment = rel(path, output_root)
    family_id = symbol
    event_no = 0
    events: List[Event] = []
    status_counter: Counter[str] = Counter()

    min_target_informative = threshold_count(args.min_target_informative, len(groups.target_tokens), default_all=True)
    min_background_informative = threshold_count(
        args.min_background_informative,
        len(groups.background_tokens),
        default_all=False,
    )
    for aln_codon_index in range(1, len(next(iter(records.values()))) // 3 + 1):
        token_states = {
            token: classify_codon(token, seq[(aln_codon_index - 1) * 3 : aln_codon_index * 3], table)
            for token, seq in records.items()
        }
        status_counter.update(s.status for s in token_states.values())
        target_stats = group_stats_for(token_states, groups.target_tokens)
        background_stats = group_stats_for(token_states, groups.background_tokens)
        outgroup_stats = group_stats_for(token_states, groups.outgroup_tokens)
        coord_state = token_states.get(groups.coordinate_reference_token)

        candidate = call_aa_variant(
            make_event_id(event_no + 1),
            family_id,
            symbol,
            source_alignment,
            aln_codon_index,
            target_stats,
            background_stats,
            outgroup_stats,
            coord_state,
            args.target_state_mode,
            args.bed_mode,
            min_target_informative,
            min_background_informative,
        )
        if candidate:
            event_no += 1
            candidate.event_id = make_event_id(event_no)
            events.append(candidate)

    aa_rows = [event_base_row(event, groups) for event in events]
    codon_rows = list(aa_rows)
    matrix_rows: List[dict] = []
    nt_rows: List[dict] = []
    for event in events:
        rows, nts = matrix_rows_for_event(event, groups, resources, output_root)
        matrix_rows.extend(rows)
        nt_rows.extend(nts)

    events_dir = output_root / "events"
    matrices_dir = output_root / "matrices"
    bed_dir = output_root / "bed"
    aa_path = events_dir / f"{symbol}.aa_events.tsv.gz"
    codon_path = events_dir / f"{symbol}.codon_events.tsv.gz"
    nt_path = events_dir / f"{symbol}.nt_changes.tsv.gz"
    matrix_path = matrices_dir / f"{symbol}.variant_matrix.tsv.gz"
    write_gzip_tsv(aa_path, AA_EVENT_COLUMNS, aa_rows)
    write_gzip_tsv(codon_path, AA_EVENT_COLUMNS, codon_rows)
    write_gzip_tsv(nt_path, NT_COLUMNS, nt_rows)
    write_gzip_tsv(matrix_path, MATRIX_COLUMNS, matrix_rows)

    target_label = sanitize_name("_".join(groups.target_tokens))
    bed_paths: List[str] = []
    total_bed_rows = 0
    for bed_class, rows in bed_rows_from_matrix(matrix_rows, symbol, target_label, args.bed_mode).items():
        bed_path = bed_dir / f"{symbol}.{target_label}.{bed_class}.bed"
        write_bed(bed_path, rows)
        bed_paths.append(rel(bed_path, output_root))
        total_bed_rows += sum(1 for _ in bed_path.open()) if bed_path.exists() else 0

    log_path = output_root / "logs" / f"{symbol}.variants.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        json.dumps(
            {
                "symbol": symbol,
                "alignment": str(path),
                "codon_status_counts": dict(status_counter),
                "target_tokens": groups.target_tokens,
                "background_tokens": groups.background_tokens,
                "outgroup_tokens": groups.outgroup_tokens,
                "exclude_tokens": groups.exclude_tokens,
                "coordinate_reference_token": groups.coordinate_reference_token,
                "coordinate_reference_role": groups.coordinate_reference_role,
                "codon_map_path": resources.codon_map_path,
                "matrix_path": resources.matrix_path,
                "events": len(events),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )

    return GeneResult(
        status="pass",
        symbol=symbol,
        family_id=family_id,
        alignment_path=rel(path, output_root),
        aa_events_path=rel(aa_path, output_root),
        codon_events_path=rel(codon_path, output_root),
        nt_changes_path=rel(nt_path, output_root),
        variant_matrix_path=rel(matrix_path, output_root),
        bed_paths=";".join(bed_paths),
        event_count=len(events),
        variant_count=sum(1 for e in events if e.event_type == "aa_variant"),
        identical_sequence_count=sum(1 for e in events if e.event_subtype == "identical_sequence"),
        divergent_sequence_count=sum(1 for e in events if e.event_subtype == "divergent_sequence"),
        substitution_count=0,
        indel_like_count=0,
        nt_change_count=len(nt_rows),
        coordinateable_rows=sum(1 for r in matrix_rows if r.get("coordinateable") == "true"),
        bed_rows=total_bed_rows,
    )


def result_row(result: GeneResult) -> dict:
    return {
        "status": result.status,
        "symbol": result.symbol,
        "family_id": result.family_id,
        "alignment_path": result.alignment_path,
        "aa_events_path": result.aa_events_path,
        "codon_events_path": result.codon_events_path,
        "nt_changes_path": result.nt_changes_path,
        "variant_matrix_path": result.variant_matrix_path,
        "bed_paths": result.bed_paths,
        "event_count": result.event_count,
        "variant_count": result.variant_count,
        "identical_sequence_count": result.identical_sequence_count,
        "divergent_sequence_count": result.divergent_sequence_count,
        "substitution_count": result.substitution_count,
        "indel_like_count": result.indel_like_count,
        "nt_change_count": result.nt_change_count,
        "coordinateable_rows": result.coordinateable_rows,
        "bed_rows": result.bed_rows,
        "reason": result.reason,
    }


def write_plain_tsv(path: Path, columns: Sequence[str], rows: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(columns), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Call target-specific amino acid variants and map coordinateable variants to BED.")
    parser.add_argument("--alignment-dir", default=str(Path.cwd() / "alignments" / "mafft_pal2nal"))
    parser.add_argument("--codon-map-dir", default="")
    parser.add_argument("--matrix-dir", default=str(Path.cwd() / "human_cds_matrices"))
    parser.add_argument("--output-dir", default=str(Path.cwd() / "variants"))
    parser.add_argument("--coordinate-reference-token", default="human")
    parser.add_argument("--reference-token", help="Deprecated alias for --coordinate-reference-token")
    parser.add_argument("--target-token", action="append", default=[])
    parser.add_argument("--target-tokens-file")
    parser.add_argument("--outgroup-token", action="append", default=[])
    parser.add_argument("--outgroup-tokens-file")
    parser.add_argument("--exclude-token", action="append", default=[])
    parser.add_argument("--exclude-tokens-file")
    parser.add_argument("--target-state-mode", choices=["uniform", "allow-diverse"], default="uniform", help="Deprecated compatibility option; v0.1.5 reports identical/divergent target states automatically")
    parser.add_argument("--min-target-non-gap", default="all", help="Deprecated alias for --min-target-informative")
    parser.add_argument("--min-target-informative", default=None, help="Minimum informative target states; valid codons and full codon gaps are informative")
    parser.add_argument("--max-target-gap-fraction-for-substitution", type=float, default=0.0)
    parser.add_argument("--min-background-non-gap", type=int, default=5, help="Deprecated alias for --min-background-informative")
    parser.add_argument("--min-background-informative", default=None, help="Minimum informative background states; valid codons and full codon gaps are informative")
    parser.add_argument("--min-background-gap-fraction-for-target-non-gap-event", type=float, default=0.8)
    parser.add_argument("--min-background-non-gap-fraction-for-target-gap-event", type=float, default=0.8)
    parser.add_argument("--min-target-non-gap-fraction", type=float, default=1.0)
    parser.add_argument("--min-target-gap-fraction", type=float, default=1.0)
    parser.add_argument("--bed-mode", choices=["auto", "all-coordinateable", "substitution-only", "none"], default="auto")
    parser.add_argument("--codon-table", choices=["universal", "vmitochondria"], default="universal")
    parser.add_argument("--symbols", nargs="*", default=[])
    parser.add_argument("--symbols-file")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--force", action="store_true")
    return parser


def normalize_args(args: argparse.Namespace) -> argparse.Namespace:
    if args.reference_token and args.coordinate_reference_token != "human" and args.reference_token != args.coordinate_reference_token:
        raise ValueError("--reference-token and --coordinate-reference-token disagree")
    if args.reference_token:
        print(
            "Warning: --reference-token is deprecated; use --coordinate-reference-token.",
            file=sys.stderr,
        )
        args.coordinate_reference_token = args.reference_token
    args.target_tokens = unique_preserve_order(split_values(args.target_token) + read_token_file(args.target_tokens_file))
    args.outgroup_tokens = unique_preserve_order(split_values(args.outgroup_token) + read_token_file(args.outgroup_tokens_file))
    args.exclude_tokens = unique_preserve_order(split_values(args.exclude_token) + read_token_file(args.exclude_tokens_file))
    args.symbols = unique_preserve_order(split_values(args.symbols) + read_token_file(args.symbols_file))
    if args.min_target_informative is None:
        args.min_target_informative = args.min_target_non_gap
    if args.min_background_informative is None:
        args.min_background_informative = str(args.min_background_non_gap)
    if args.bed_mode == "auto":
        args.bed_mode = "all-coordinateable"
    elif args.bed_mode == "substitution-only":
        print(
            "Warning: --bed-mode substitution-only is deprecated in v0.1.5; writing coordinateable variant rows.",
            file=sys.stderr,
        )
        args.bed_mode = "all-coordinateable"
    if not args.codon_map_dir:
        aln_dir = Path(args.alignment_dir)
        if (aln_dir / "maps").is_dir():
            args.codon_map_dir = str(aln_dir / "maps")
        elif (aln_dir.parent / "maps").is_dir():
            args.codon_map_dir = str(aln_dir.parent / "maps")
        else:
            args.codon_map_dir = str(Path.cwd() / "alignments" / "maps")
    return args


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_parser()
    args = normalize_args(parser.parse_args(argv))
    start = time.time()
    output_root = Path(args.output_dir).resolve()
    if output_root.exists() and any(output_root.iterdir()) and not args.force:
        raise FileExistsError(f"{output_root} is not empty. Use --force to overwrite/reuse outputs.")
    for subdir in ["events", "matrices", "bed", "logs"]:
        (output_root / subdir).mkdir(parents=True, exist_ok=True)
    alignment_files = collect_alignment_files(Path(args.alignment_dir).resolve(), args.symbols, args.limit)
    preflight = {
        "alignment_dir": str(Path(args.alignment_dir).resolve()),
        "codon_map_dir": str(Path(args.codon_map_dir).resolve()),
        "matrix_dir": str(Path(args.matrix_dir).resolve()),
        "output_dir": str(output_root),
        "coordinate_reference_token": args.coordinate_reference_token,
        "target_tokens": args.target_tokens,
        "outgroup_tokens": args.outgroup_tokens,
        "exclude_tokens": args.exclude_tokens,
        "target_state_mode": args.target_state_mode,
        "min_target_informative": args.min_target_informative,
        "min_background_informative": args.min_background_informative,
        "bed_mode": args.bed_mode,
        "selected_alignment_files": [str(p) for p in alignment_files],
    }
    (output_root / "preflight.json").write_text(json.dumps(preflight, indent=2, sort_keys=True) + "\n")

    results: List[GeneResult] = []
    failed: List[GeneResult] = []
    for path in alignment_files:
        symbol = symbol_from_alignment(path)
        log(f"Calling variants for {symbol}")
        try:
            results.append(process_gene(path, args, output_root))
        except Exception as exc:
            failed.append(
                GeneResult(
                    status="fail",
                    symbol=symbol,
                    family_id=symbol,
                    alignment_path=str(path),
                    reason=f"{type(exc).__name__}: {exc}",
                )
            )
    columns = list(result_row(GeneResult("pass", "", "", "")).keys())
    write_plain_tsv(output_root / "manifest.tsv", columns, [result_row(r) for r in results])
    write_plain_tsv(output_root / "failed.tsv", columns, [result_row(r) for r in failed])
    merged_bed_path, merged_bed_rows = write_merged_bed(output_root, results)
    summary = {
        "status": "pass" if not failed else "partial",
        "version_stage": "v0.1.5",
        "families_processed": len(results),
        "families_failed": len(failed),
        "variant_calling_mode": "aa_mutual_exclusive",
        "target_state_mode": args.target_state_mode,
        "min_target_informative": args.min_target_informative,
        "min_background_informative": args.min_background_informative,
        "bed_mode": args.bed_mode,
        "coordinate_reference_token": args.coordinate_reference_token,
        "target_tokens": args.target_tokens,
        "outgroup_tokens": args.outgroup_tokens,
        "exclude_tokens": args.exclude_tokens,
        "event_count": sum(r.event_count for r in results),
        "variant_count": sum(r.variant_count for r in results),
        "identical_sequence_count": sum(r.identical_sequence_count for r in results),
        "divergent_sequence_count": sum(r.divergent_sequence_count for r in results),
        "substitution_count": 0,
        "indel_like_count": 0,
        "nt_change_count": sum(r.nt_change_count for r in results),
        "coordinateable_rows": sum(r.coordinateable_rows for r in results),
        "bed_rows": sum(r.bed_rows for r in results),
        "merged_bed_path": merged_bed_path,
        "merged_bed_rows": merged_bed_rows,
        "runtime_seconds": round(time.time() - start, 3),
    }
    (output_root / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
