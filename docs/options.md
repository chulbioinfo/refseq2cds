# Command Options

This page summarizes the options most users need after installation.

## `refseq2cds run`

Build CDS FASTA files with either strict all-species singleton parsing or a
single reference-gene present-species query.

```bash
refseq2cds run --steps all --with-matrices
```

Important options:

| Option | Meaning |
|---|---|
| `--manifest PATH` | Species/assembly table for the analysis. Defaults to `./config/species_manifest.tsv` if present, otherwise the packaged default manifest. |
| `--orthology-mode strict_singleton` | Default mode. Keep only families with exactly one GeneID and one valid CDS in every manifest species. |
| `--orthology-mode reference_gene_1to1_present_species` | Query one reference gene and keep only manifest species with an unambiguous 1:1 orthology relationship to that gene. |
| `--reference-taxid TAXID` | Species whose gene symbols name FASTA files. Defaults to human `9606`. |
| `--reference-symbol SYMBOL` | Reference gene symbol for `reference_gene_1to1_present_species`, for example `FOXP2`. |
| `--reference-gene-id GENEID` | Reference NCBI GeneID for `reference_gene_1to1_present_species`. Use this when a symbol is ambiguous. |
| `--min-sequences N` | Minimum output sequence count for reference-gene mode. Defaults to `2`. |
| `--exclude-reference` | Omit the reference species sequence from reference-gene mode outputs. Usually leave this off. |
| `--with-matrices` | After FASTA generation, build human CDS-to-genome coordinate matrices if human is in the manifest. |
| `--input-root PATH` | Directory containing pre-existing `ncbi_bulk/` and `assembly_packages/`. Mainly for offline tests or managed HPC downloads. |
| `--output-root PATH` | Directory where generated indexes, FASTA files, reports, and matrices are written. |
| `--offline` | Use local fixtures under `--input-root` and do not download from NCBI. |
| `--force` | Rebuild stage outputs even if files already exist. |

Reference-gene example:

```bash
refseq2cds run \
  --steps all \
  --orthology-mode reference_gene_1to1_present_species \
  --reference-taxid 9606 \
  --reference-symbol FOXP2 \
  --with-matrices \
  --force
```

This mode excludes species with no NCBI Gene orthology edge, one-reference to
many-target ambiguity, many-reference to one-target ambiguity, missing locked
GCF annotation, non-protein-coding or mitochondrial genes, missing CDS, or CDS
QC failure. The output species set can therefore be smaller than the manifest.
For large manifests, the pipeline writes
`reports/reference_gene_download_scope.json` before assembly package downloads;
use it to estimate how many GCF annotation packages the query will require.

## `refseq2cds test`

Run an offline smoke test.

```bash
refseq2cds test --example mini
```

The mini example checks orthology graph filtering, CDS QC, FASTA writing, and
human coordinate matrix generation without downloading data from NCBI.

## `refseq2cds verify`

Check generated outputs.

```bash
refseq2cds verify --output-root /path/to/run --manifest /path/to/species_manifest.tsv
```

If `--output-root` is omitted, the current directory is used.

## `refseq2cds manifest-from-gcf`

Create a manifest from RefSeq `GCF_*` accessions.

```bash
refseq2cds manifest-from-gcf --inputfile gcfs.txt --download-tools --force
```

The generated manifest records each species token, taxid, scientific name,
common name, GCF accession, and optional outgroup flag. Inputs such as
`GCF_037993035.2-RS_2025_03` or `GCF_964374335.1/` are normalized to the base
RefSeq assembly accession before querying NCBI.

## `refseq2cds build-matrices`

Rebuild human coordinate matrices from existing FASTA and selection outputs.

```bash
refseq2cds build-matrices --force
```

## `refseq2cds align`

Build lightweight tree-free codon alignments with MAFFT protein MSA and PAL2NAL
back-translation.

```bash
refseq2cds align --mode mafft-pal2nal --jobs 4 --threads-per-mafft 2
```

This is a baseline alignment mode. It does not replace tree-guided codon
alignment tools such as PAGAN2 or PRANK.

## `refseq2cds variants`

Call target-group-specific comparative coding events from codon-aware
alignments and map coordinateable events to BED.

```bash
refseq2cds variants \
  --alignment-dir alignments/mafft_pal2nal \
  --codon-map-dir alignments/mafft_pal2nal/maps \
  --matrix-dir human_cds_matrices \
  --coordinate-reference-token human \
  --target-token human \
  --min-background-non-gap 5 \
  --force
```

Important options:

| Option | Meaning |
|---|---|
| `--coordinate-reference-token TOKEN` | Species token used only for CDS/genome coordinate mapping. It is not a privileged event-calling comparator. Defaults to `human`. |
| `--target-token TOKEN` | Target species or group member. May be supplied more than once. |
| `--target-tokens-file PATH` | File with target tokens. |
| `--outgroup-token TOKEN` | Biological outgroup token excluded from target-vs-background event calling. May be supplied more than once. |
| `--exclude-token TOKEN` | Token excluded from event calling for quality or analysis reasons. |
| `--target-state-mode uniform` | Default. Target tokens must share the same valid amino acid state for substitution calls. |
| `--target-state-mode allow-diverse` | Allows multiple target states if the target state set and background state set are disjoint. |
| `--bed-mode auto` | Write all coordinateable event classes to BED. Alias of `all-coordinateable`. |
| `--bed-mode substitution-only` | Write only nonsynonymous substitution-derived BED rows. |
| `--bed-mode none` | Write matrix/event tables only. |

The `variants` command compares target states against background states after
removing outgroup and excluded tokens. The coordinate reference species is used
only after an event has already been called, when alignment positions are mapped
through the coordinate-reference codon map and CDS-to-genome matrix.
