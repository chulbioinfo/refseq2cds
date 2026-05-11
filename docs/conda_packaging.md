# Conda Packaging Notes

`refseq2cds` includes a draft conda recipe in `conda-recipe/meta.yaml`.
This is intended to make a future Bioconda submission easier.

## Local conda-build test

Install conda-build tools:

```bash
mamba install -c conda-forge -c bioconda conda-build conda-verify
```

Build locally:

```bash
conda build conda-recipe
```

For local builds, use a clean clone or remove large generated directories first.
The committed source distribution excludes `raw/`, `fastas/`, indexes, reports
and other generated outputs, but a local `source: path: ..` conda build may
still scan the working tree before packaging.

The recipe test runs:

```bash
refseq2cds --help
refseq2cds test --example mini
```

The mini test is offline and uses packaged fixtures, so it should not contact
NCBI during conda-build.

## Bioconda submission checklist

Before opening a Bioconda pull request:

- create a GitHub release tag, for example `v0.1.2`;
- upload a source archive or use the GitHub release tarball URL;
- replace `source: path: ..` with a release `url` and `sha256`;
- confirm all run dependencies exist on `conda-forge` or `bioconda`;
- confirm `refseq2cds test --example mini` passes inside the conda-build test
  environment;
- keep command-line dependencies (`ncbi-datasets-cli`, `mafft`, `pal2nal`) as
  conda run requirements because they are used by normal workflows.

## Why package data is included

The CLI dispatches to scripts under `workflow/scripts` and the reviewer mini
fixture under `examples/mini`. Non-editable wheel and conda installs therefore
need those files under `share/refseq2cds`. The CLI searches the repository first
and then the installed share directory.
