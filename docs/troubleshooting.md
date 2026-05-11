# Troubleshooting

## `datasets` is missing

Run:

```bash
refseq2cds download-tools
```

or install NCBI Datasets CLI with conda:

```bash
mamba install -c conda-forge ncbi-datasets-cli
```

Then pass the executable explicitly if needed:

```bash
refseq2cds run --datasets /path/to/datasets
```

## Windows support

Windows is not currently supported. Use macOS, Linux, WSL with caution, or a
Linux container/HPC environment. The automatic NCBI Datasets downloader supports
macOS and Linux only.

## Changing the manifest

Changing `config/species_manifest.tsv` changes the biological dataset. Rebuild
with:

```bash
refseq2cds run --steps all --force --with-matrices
```

For separate experiments, prefer `--output-root` so old and new results do not
mix.

## Human matrices are skipped

Human matrices require taxid `9606` in the manifest and `--with-matrices`:

```bash
refseq2cds run --steps all --with-matrices
```

If the manifest uses a non-human reference species, human matrices can still be
created as long as human is included in the species set.

## MAFFT or PAL2NAL is missing

Install the optional alignment dependencies:

```bash
mamba install -c bioconda -c conda-forge mafft pal2nal
```

or pass explicit paths:

```bash
refseq2cds align --mafft /path/to/mafft --pal2nal /path/to/pal2nal.pl
```

## Verify the installation without downloading NCBI data

Run:

```bash
refseq2cds test --example mini
```

This should finish in seconds and produce two FASTA files plus two human
coordinate matrices in a temporary output directory.
