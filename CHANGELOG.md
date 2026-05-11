# Changelog

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
