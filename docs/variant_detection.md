# Variant Detection and Coordinate Mapping

`refseq2cds variants` runs after codon-aware alignment. It detects
target-group-specific comparative coding events by comparing target states
against background states.

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

Default mode:

```bash
--target-state-mode uniform
```

This requires target tokens to share one valid amino acid state. The
`allow-diverse` mode permits multiple target states if the target state set and
background state set are disjoint.

## Outputs

```text
variants/events/{SYMBOL}.aa_events.tsv.gz
variants/events/{SYMBOL}.codon_events.tsv.gz
variants/events/{SYMBOL}.nt_changes.tsv.gz
variants/matrices/{SYMBOL}.variant_matrix.tsv.gz
variants/bed/{SYMBOL}.{TARGET_SET}.target_exclusive_substitutions.bed
variants/bed/{SYMBOL}.{TARGET_SET}.target_non_gap_background_gap.bed
variants/bed/{SYMBOL}.{TARGET_SET}.target_gap_background_non_gap.bed
variants/summary.json
```

BED files contain only events where the coordinate reference has a real
mappable base. Events without a coordinate-reference base remain in the variant
matrix with a `coordinate_status` value explaining why no BED row was written.

## BED Interpretation for Indel-like Events

The two indel-like classes are alignment-relative, so BED conversion depends on
which species is used as the coordinate reference:

- `target_non_gap_background_gap`: the target has a codon where the background
  is mostly gapped. This can be written to BED when the coordinate reference has
  that target/base, for example when human is both `--target-token human` and
  `--coordinate-reference-token human`.
- `target_gap_background_non_gap`: the target is gapped where the background has
  a codon. This can be written to BED when the coordinate reference is in the
  background and therefore has the mappable base.
- if the coordinate reference is gapped at the event, no BED interval is written
  because there is no reference-genome base to anchor.

## Limitations

v0.1.4 does not infer ancestral states, true insertion/deletion polarity, or
derived versus ancestral status. Indel-like event names are alignment-relative:
`target_non_gap_background_gap` and `target_gap_background_non_gap`.
