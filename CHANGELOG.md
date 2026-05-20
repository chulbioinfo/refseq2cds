# Changelog

## v0.1.5 - 2026-05-20

### Added

- `refseq2cds --version` for quick installed-version checks.
- Optional `refseq2cds verify --alignment-dir ...` and `--variant-dir ...`
  checks so end-to-end runs can validate FASTA/matrix outputs together with
  MAFFT+PAL2NAL alignment and variant-call directories. Relative paths are
  resolved under `--output-root`.
- `refseq2cds verify --matrix-rows none` now supports FASTA/report-only runs
  that were generated without `--with-matrices`.

### Changed

- Project version bumped to `0.1.5`.
- `refseq2cds variants` now defaults to `alignments/mafft_pal2nal`, matching
  the default output directory from `refseq2cds align`.
- README and option documentation were reorganized around the practical
  `manifest -> run -> verify -> align -> variants -> verify` analysis flow.
- Reworked `refseq2cds variants` to call codon-site variants only when target
  and background amino acid state sets are mutually exclusive.
- Full codon gaps now participate as amino acid state `-`; gap fraction
  thresholds no longer create separate insertion/deletion calls.
- Uniform target states are reported as `identical_sequence`; diverse target
  states are reported as `divergent_sequence`.
- Coordinateable BED output is unified under `bed_event_class=variant` and
  `{SYMBOL}.{TARGET_SET}.variant.bed`.
- Added `variants/merged.bed`, a coordinate-sorted concatenation of all
  per-gene variant BED rows for downstream interval analyses.

## v0.1.4 - 2026-05-11

### Added

- `refseq2cds variants` for target-group versus background-group comparative coding event detection from codon-aware alignments.
- Target-specific amino acid substitution and alignment-relative indel-like event outputs.
- Nonsynonymous nucleotide component extraction for target-exclusive amino acid substitution events.
- Coordinate-reference mapping from alignment codon positions through codon maps and CDS-to-genome matrices to BED6.
- `--target-state-mode uniform` default and `--target-state-mode allow-diverse` for diverse target groups.
- Outgroup and exclude token handling for event calling.
- Variant matrix, per-gene event tables, BED outputs, summary, failure log, and tests.

### Changed

- Project version bumped to `0.1.4`.
- The coordinate reference token is explicitly used only for coordinate mapping and is not a privileged biological comparator.

## v0.1.3 - 2026-05-11

### Added

- `reference_gene_1to1_present_species` orthology mode for single reference-gene queries such as human `FOXP2`.
- Reference-gene graph filtering that keeps species with one unambiguous target GeneID and excludes no-edge, 1:M, M:1, M:M, non-coding, mitochondrial, missing-annotation, and invalid-CDS cases.
- Per-query output directories under `results/reference_gene_1to1/{REFERENCE_SYMBOL}/` with candidate, pass, rejection, CDS selection, QC, FASTA, metadata, and summary files.
- CLI options for `--reference-symbol`, `--reference-gene-id`, `--min-sequences`, and `--exclude-reference`.
- Tests for present-species reference-gene extraction, missing target species, and one-reference-to-many-target rejection.

### Changed

- Project version bumped to `0.1.3`.
- FASTA verification now accepts subset species only for reference-gene present-species outputs while keeping strict full-manifest checks for strict singleton mode.
- Run configuration snapshots now prevent accidental reuse of old outputs when orthology mode or reference-gene settings change without `--force`.
- `pyproject.toml` now uses SPDX-style `license = "MIT"` for cleaner package builds.

## v0.1.2 - 2026-05-11

### Added

- Offline `refseq2cds test --example mini` smoke test covering strict singleton filtering, CDS QC, FASTA writing, and human coordinate matrix generation.
- Packaged `examples/mini/` fixture for reviewer/developer verification without NCBI downloads.
- Pytest tests and GitHub Actions CI workflow.
- `environment.yml` for reproducible conda/mamba setup.
- Draft `conda-recipe/meta.yaml` for future conda/Bioconda packaging.
- `CITATION.cff` and packaging/troubleshooting documentation.

### Changed

- Project version bumped to `0.1.2`.
- Pipeline drivers now support `--manifest`, `--input-root`, `--output-root`, and offline fixture execution.
- Installable package metadata now includes workflow scripts and mini fixtures under `share/refseq2cds`.

## v0.1.1 - 2026-05-11

### Added

- `refseq2cds align --mode mafft-pal2nal` for lightweight tree-free codon alignment.
- MAFFT protein-alignment generation from each per-family CDS FASTA.
- PAL2NAL codon back-translation from the MAFFT protein MSA and original CDS FASTA.
- Alignment manifest, failure log, summary JSON, per-gene logs, and optional token-specific codon maps.
- README documentation for optional MAFFT/PAL2NAL dependencies and alignment usage.

### Changed

- Project version bumped to `0.1.1`.
- README now distinguishes tree-free MAFFT+PAL2NAL baseline alignments from tree-guided PAGAN2/PRANK analyses.

## v0.1.0 - 2026-05-11

Initial release of `refseq2cds`.

### Added

- Assembly-exact RefSeq CDS extraction from user-provided `GCF_*` accessions.
- Automatic manifest generation from GCF accessions with NCBI taxid and species metadata.
- Strict N-way singleton ortholog detection from NCBI Gene FTP `gene_orthologs.gz`.
- Protein-coding, non-mitochondrial filtering.
- Representative CDS selection with RefSeq transcript/protein metadata.
- Per-family `{REFERENCE_SYMBOL}.fasta` and metadata output.
- Optional human CDS nucleotide to genomic coordinate matrices for downstream BED/UCSC workflows.
- Configurable reference species for output symbols and FASTA filenames.
- Verification and summary commands.
- Beginner-friendly README covering GCF accessions, taxids, manifests, reference species, and supported platforms.

### Notes

- Supported platforms: macOS and Linux.
- Windows is not currently supported.
- Intended biological scope: vertebrate RefSeq workflows only.
- Alignment is not performed by this package yet; FASTA outputs are ready for PAGAN2/PRANK or other codon-aware aligners.
