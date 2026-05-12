# Conda Packaging Notes

`refseq2cds` includes a Bioconda-ready recipe in `conda-recipe/meta.yaml`.
The recipe should always be pinned to a GitHub release tarball and SHA256. The
checked-in recipe currently targets the already tagged v0.1.2 release:

```text
https://github.com/chulbioinfo/refseq2cds/archive/refs/tags/v0.1.2.tar.gz
sha256: 945b7876eb314c5809357bb4067d3934a706641b9d24987b2026a6d7600dd474
```

## Local conda-build test

Install conda-build tools:

```bash
mamba install -c conda-forge -c bioconda conda-build
```

Build locally:

```bash
conda build conda-recipe
```

The recipe now builds from the release tarball, not from the local working tree,
so ignored generated directories such as `raw/`, `fastas/`, indexes, reports,
and matrices are not part of the conda source.

The recipe test runs:

```bash
refseq2cds --help
refseq2cds test --example mini
```

The mini test is offline and uses packaged fixtures, so it should not contact
NCBI during conda-build.

## Bioconda submission checklist

Before opening a Bioconda pull request, verify:

- the GitHub tag and release for the target version exist;
- `conda-recipe/meta.yaml` has the release tarball URL and SHA256;
- all run dependencies exist on `conda-forge` or `bioconda`;
- confirm `refseq2cds test --example mini` passes inside the conda-build test
  environment;
- keep command-line dependencies (`ncbi-datasets-cli`, `mafft`, `pal2nal`) as
  conda run requirements because they are used by normal workflows.

For a new release such as v0.1.4, tag and push the GitHub release first, then
download `https://github.com/chulbioinfo/refseq2cds/archive/refs/tags/vX.Y.Z.tar.gz`,
compute the SHA256, update `conda-recipe/meta.yaml`, and run the lint/build
checks before opening or updating the Bioconda PR.

Dependency channels checked:

| Package | Channel |
|---|---|
| `ncbi-datasets-cli` | `conda-forge` |
| `mafft` | `bioconda` |
| `pal2nal` | `bioconda` |
| Python libraries (`biopython`, `dendropy`, `jinja2`, `networkx`, `pandas`, `pyarrow`) | `conda-forge`/`bioconda` |

## Bioconda PR workflow

Fork and clone the Bioconda recipes repository:

```bash
git clone https://github.com/<YOUR_GITHUB_USERNAME>/bioconda-recipes.git
cd bioconda-recipes
git remote add upstream https://github.com/bioconda/bioconda-recipes.git
git checkout -b add-refseq2cds
```

Copy the recipe:

```bash
mkdir -p recipes/refseq2cds
cp /Users/openhl/Documents/NCBI_ortholog/conda-recipe/meta.yaml recipes/refseq2cds/meta.yaml
```

Run local lint/build if you have Bioconda tooling installed:

```bash
bioconda-utils lint recipes config.yml --packages refseq2cds
bioconda-utils build recipes config.yml --packages refseq2cds
```

Then commit, push, and open a pull request:

```bash
git add recipes/refseq2cds/meta.yaml
git commit -m "Add refseq2cds"
git push origin add-refseq2cds
```

## Why package data is included

The CLI dispatches to scripts under `workflow/scripts` and the reviewer mini
fixture under `examples/mini`. Non-editable wheel and conda installs therefore
need those files under `share/refseq2cds`. The CLI searches the repository first
and then the installed share directory.
