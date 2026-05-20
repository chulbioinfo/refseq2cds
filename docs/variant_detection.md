# Variant Detection and Coordinate Mapping

`refseq2cds variants` runs after codon-aware alignment. In v0.1.5 it detects
target-group-specific amino acid variants by comparing target amino acid states
against background amino acid states at each alignment codon.

The coordinate reference species is used only for coordinate mapping:

```text
alignment codon index
→ coordinate-reference codon map
→ coordinate-reference CDS nucleotide position
→ CDS-to-genome matrix
→ BED6 coordinate
```

It is not a privileged biological comparator during event calling.

## Event Calling

For each alignment codon site:

```text
background = all tokens - target tokens - outgroup tokens - excluded tokens
```

Outgroup tokens are reported but excluded from the decision.

Full codon gaps (`---`) are treated as amino acid state `-`. Valid codons and
full codon gaps are informative. Partial gaps, ambiguous codons, invalid
codons, and stop codons are counted for QC but do not create v0.1.5 calls by
default.

A site is called only when target and background states are mutually exclusive:

```text
target_state_set is not empty
background_state_set is not empty
target_state_set ∩ background_state_set == ∅
```

Target states define the subtype:

```text
len(target_state_set) == 1  → identical_sequence
len(target_state_set) >= 2  → divergent_sequence
```

The caller does not label primary outputs as substitution, insertion, or
deletion. All coordinateable calls are emitted with `bed_event_class=variant`.

## Outputs

```text
variants/events/{SYMBOL}.aa_events.tsv.gz
variants/events/{SYMBOL}.codon_events.tsv.gz
variants/events/{SYMBOL}.nt_changes.tsv.gz
variants/matrices/{SYMBOL}.variant_matrix.tsv.gz
variants/bed/{SYMBOL}.{TARGET_SET}.variant.bed
variants/merged.bed
variants/summary.json
```

BED files contain only events where the coordinate reference has a real
mappable base. Events without a coordinate-reference base remain in the variant
matrix with a `coordinate_status` value explaining why no BED row was written.
`variants/merged.bed` concatenates all per-gene `*.variant.bed` files and sorts
rows by chromosome, start, end, and name. `summary.json` records both
`bed_rows` and `merged_bed_rows`; these counts should match.

## Coordinate Mapping

After a variant is called, genomic coordinate conversion is attempted only when
the coordinate reference has a valid non-gap codon at the same alignment site.
If the coordinate reference is gapped, partially gapped, ambiguous, missing, or
not present in the codon map/matrix, the event remains matrix-only.

Common `coordinate_status` values include:

```text
mapped_coordinate_reference_base
not_mapped_coordinate_reference_gap
not_mapped_coordinate_reference_partial_gap
not_mapped_coordinate_reference_ambiguous_nt
not_mapped_coordinate_reference_stop_codon
not_mapped_missing_codon_map
not_mapped_missing_reference_matrix
not_mapped_no_genomic_coordinate
```

## Limitations

v0.1.5 does not infer ancestral states, true insertion/deletion polarity, derived
status, or ancestral direction. It deliberately keeps event calling at the
amino acid state level and uses genomic mapping only after the mutually
exclusive amino acid variant has already been detected.
