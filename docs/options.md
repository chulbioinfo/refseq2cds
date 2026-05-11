# Command Options

This page summarizes the options most users need after installation.

## `refseq2cds run`

Build strict singleton CDS FASTA files.

```bash
refseq2cds run --steps all --with-matrices
```

Important options:

| Option | Meaning |
|---|---|
| `--manifest PATH` | Species/assembly table for the analysis. Defaults to `./config/species_manifest.tsv` if present, otherwise the packaged default manifest. |
| `--reference-taxid TAXID` | Species whose gene symbols name FASTA files. Defaults to human `9606`. |
| `--with-matrices` | After FASTA generation, build human CDS-to-genome coordinate matrices if human is in the manifest. |
| `--input-root PATH` | Directory containing pre-existing `ncbi_bulk/` and `assembly_packages/`. Mainly for offline tests or managed HPC downloads. |
| `--output-root PATH` | Directory where generated indexes, FASTA files, reports, and matrices are written. |
| `--offline` | Use local fixtures under `--input-root` and do not download from NCBI. |
| `--force` | Rebuild stage outputs even if files already exist. |

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
common name, GCF accession, and optional outgroup flag.

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
