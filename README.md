# refseq2cds

Assembly-exact RefSeq ortholog CDS FASTA, alignment, and variant workflow.

Version: `0.1.5`

Author: Chul Lee (chul.bioinfo@gmail.com)

`refseq2cds` is an assembly-exact NCBI RefSeq ortholog CDS FASTA builder with
optional coordinate matrices, codon alignments, and target-specific variant
outputs.

This repository builds assembly-exact ortholog CDS FASTA files from NCBI RefSeq
assembly annotations and NCBI Gene FTP orthology edges. It supports two
orthology parsing modes: strict N-way singleton extraction across every locked
species, and a reference-gene present-species mode for queries such as human
`FOXP2` where only species with an unambiguous 1:1 relationship to the query
gene are retained. It also creates human CDS-position to genomic-position
matrices so codon/site-level results can later be converted to UCSC Genome
Browser BED coordinates. It can also create lightweight tree-free codon
alignments with a MAFFT protein MSA followed by PAL2NAL back-translation.
After codon alignment, `refseq2cds variants` can detect target-group-specific
amino acid variants and write coordinateable `variant` events as BED.

`refseq2cds` is not a de novo orthology inference program and it is not just a
CDS downloader. It treats NCBI Gene orthology records as input evidence, applies
auditable graph filters, selects CDS records from exact RefSeq `GCF_*` assembly
annotation packages, records rejection reasons, and emits FASTA, metadata,
report, and coordinate-matrix outputs.

## Scope

This workflow is intended for vertebrates only. Although NCBI also distributes
some orthology data outside vertebrates, this package is built and tested for
vertebrate RefSeq genome annotation packages, standard nuclear coding genes, and
vertebrate-style comparative CDS analyses. Do not use it as-is for insects,
plants, fungi, bacteria, or other non-vertebrate projects.

The repository ships with a default 14 primate/Dermoptera manifest, but you can
also provide your own list of `GCF_*` RefSeq assembly accessions. Human
taxid `9606` is required only if you want human CDS-to-genome coordinate
matrices.

## Key Ideas For First-Time Users

### What is a GCF accession?

A `GCF_*` accession is an NCBI RefSeq genome assembly accession. In this
pipeline it answers the question: "Which exact RefSeq assembly annotation should
the CDS sequences come from?"

Example:

```text
GCF_009914755.1  # human T2T-CHM13v2.0 RefSeq assembly
GCF_028858775.2  # chimpanzee RefSeq assembly
```

This pipeline uses `datasets download genome accession <GCF>` so the CDS, RNA,
protein, and GFF3 annotation all come from that exact assembly package.

### What is a taxid?

A taxid is the numeric NCBI Taxonomy identifier for a species or taxon. For
example:

```text
9606  Homo sapiens
9598  Pan troglodytes
```

NCBI Gene orthology is keyed by `taxid` and `GeneID`, so this pipeline uses
taxids internally to decide whether an ortholog component contains exactly one
gene from each species.

You can look up a taxid from a species name with NCBI Datasets:

```bash
refseq2cds download-tools
bin/datasets summary taxonomy taxon "Homo sapiens" \
  --report ids_only \
  --as-json-lines
```

Expected output:

```json
{"query":["Homo sapiens"],"taxonomy":{"tax_id":9606}}
```

If you start from `GCF_*` assembly accessions, you usually do not need to look
up taxids manually. `refseq2cds manifest-from-gcf` queries NCBI and fills them
in for you.

### What is the manifest?

The manifest is the fixed table of species and assemblies for one analysis. It
is stored at:

```text
config/species_manifest.tsv
```

It tells the pipeline:

- which species are included
- which NCBI taxid belongs to each species
- which exact RefSeq `GCF_*` assembly to download
- what short token should be used in FASTA headers
- optionally, which species is an outgroup for later tree/alignment work

Example:

```text
token   taxid  scientific_name  common_name  gcf_accession    outgroup
human   9606   Homo sapiens     human        GCF_009914755.1  false
chimpanzee  9598  Pan troglodytes  chimpanzee  GCF_028858775.2  false
```

You usually do not need to write this by hand. Put your GCF accessions in a file
and run:

```bash
refseq2cds manifest-from-gcf --inputfile gcfs.txt --force
```

That command queries NCBI and fills in taxids and species names automatically.

### What is the reference species?

The reference species is the species whose gene symbols are used to name
ortholog families and FASTA files. By default this is human:

```bash
refseq2cds run --steps all --reference-taxid 9606
```

If you want chimpanzee symbols and filenames instead:

```bash
refseq2cds run --steps all --reference-taxid 9598
```

This does not change the orthology graph itself. It changes which species must
be present as the naming anchor and which species provides the output symbols.
Human CDS genomic matrices are a separate optional output and require human
taxid `9606` to be present in the manifest.

### What is a strict singleton ortholog?

In this tool, a strict singleton family means:

- the NCBI orthology graph contains exactly one gene from each manifest taxid
- there are no missing species
- there are no duplicated genes/paralogs within any species for that component
- all retained genes are protein-coding and non-mitochondrial
- every species has one valid CDS selected from its frozen RefSeq assembly

This conservative filter is why some genes are rejected.

### What is reference-gene present-species mode?

Strict singleton mode is best when you want a genome-wide matrix where every
output FASTA has the same species set. Sometimes you instead care about one
reference gene, such as human `FOXP2`, and want every manifest species that has
a clean ortholog to that gene.

Use:

```bash
refseq2cds run \
  --steps all \
  --orthology-mode reference_gene_1to1_present_species \
  --reference-taxid 9606 \
  --reference-symbol FOXP2 \
  --with-matrices
```

In this mode, each non-reference species is evaluated independently:

- if the species has no NCBI Gene orthology edge to the reference gene, it is
  excluded
- if the reference gene maps to multiple genes in that species, the species is
  excluded
- if the target gene maps back to multiple reference-species genes, the species
  is excluded
- if the target gene is not present in the locked `GCF_*` annotation package,
  is non-protein-coding, mitochondrial, lacks a valid CDS, or fails CDS QC, it
  is excluded

The resulting FASTA therefore contains the reference sequence plus only the
species that passed the unambiguous 1:1 graph and CDS checks. It does not force
all manifest species to be present.

For large manifests, reference-gene mode first uses the NCBI Gene bulk files to
estimate the graph-level 1:1 target species and writes
`reports/reference_gene_download_scope.json`. Only the reference species and
graph-passing target species need assembly packages. This saves time for sparse
queries, but broadly conserved genes such as `FOXP2` can still require hundreds
of large RefSeq annotation packages.

## What It Produces

Main CDS outputs:

```text
fastas/{REFERENCE_SYMBOL}.fasta
fastas/{REFERENCE_SYMBOL}.meta.tsv
fastas/manifest.tsv
```

Reference-gene present-species mode additionally writes a per-query directory:

```text
results/reference_gene_1to1/{REFERENCE_SYMBOL}/
├── reference_gene.tsv
├── ortholog_candidates.tsv
├── ortholog_pass_graph_1to1.tsv
├── ortholog_rejected_graph.tsv
├── ortholog_coverage.json
├── selected_cds.tsv
├── cds_qc.tsv
├── cds_pass.tsv
├── cds_rejected.tsv
├── {REFERENCE_SYMBOL}.reference_1to1.cds.fasta
├── {REFERENCE_SYMBOL}.reference_1to1.meta.tsv
├── {REFERENCE_SYMBOL}.reference_1to1.rejected.tsv
├── summary.json
└── summary.html
```

Optional human coordinate matrix outputs, when `--with-matrices` is used and
human taxid `9606` is present in the manifest:

```text
human_cds_matrices/{HUMAN_SYMBOL}.human_cds_genomic_matrix.tsv.gz
human_cds_matrices/manifest.tsv
human_cds_matrices/failed.tsv
```

Summary reports:

```text
reports/summary.json
reports/summary.html
reports/human_cds_position_matrices.summary.json
```

Optional MAFFT+PAL2NAL codon alignment outputs:

```text
alignments/mafft_pal2nal/codon/{REFERENCE_SYMBOL}.codon.fasta
alignments/mafft_pal2nal/aa/{REFERENCE_SYMBOL}.aa.fasta
alignments/mafft_pal2nal/aa_aligned/{REFERENCE_SYMBOL}.aa.aln.fasta
alignments/mafft_pal2nal/manifest.tsv
alignments/mafft_pal2nal/failed.tsv
alignments/mafft_pal2nal/summary.json
```

Optional target-specific variant outputs:

```text
variants/events/{REFERENCE_SYMBOL}.aa_events.tsv.gz
variants/events/{REFERENCE_SYMBOL}.codon_events.tsv.gz
variants/events/{REFERENCE_SYMBOL}.nt_changes.tsv.gz
variants/matrices/{REFERENCE_SYMBOL}.variant_matrix.tsv.gz
variants/bed/{REFERENCE_SYMBOL}.{TARGET_SET}.variant.bed
variants/merged.bed
variants/manifest.tsv
variants/summary.json
```

In strict singleton mode, each FASTA contains exactly one CDS sequence per
locked species. In reference-gene present-species mode, each FASTA contains the
reference sequence plus the subset of manifest species that pass the
unambiguous 1:1 query-gene checks. Headers are the locked species tokens, such
as:

```text
>human
>chimpanzee
>bonobo
...
>Philippine_flying_lemur
```

## Pipeline Logic

1. Freeze NCBI Gene FTP files:
   - `gene_orthologs.gz`
   - `gene_info.gz`
2. Download exact annotation packages for the frozen `GCF_*` assemblies:
   - CDS FASTA
   - protein FASTA
   - RNA FASTA
   - GFF3
   - sequence report
3. Parse assembly-exact CDS/transcript/protein/GeneID indexes.
4. Build the locked-species pairwise orthology graph from `gene_orthologs.gz`.
5. Apply one orthology parsing mode:
   - strict singleton: keep only connected components with exactly one GeneID
     per locked taxid
   - reference-gene present-species: keep species with exactly one target GeneID
     linked to the query reference GeneID and no reverse M:1 ambiguity
6. Reject MT and non-protein-coding genes.
7. Select one assembly-derived representative CDS per gene/species.
8. QC CDS sequences:
   - length multiple of 3
   - no internal stop codons
   - valid DNA alphabet
   - terminal stop removed as CDS normalization
9. Write one `{REFERENCE_SYMBOL}.fasta` per accepted singleton family.
10. Optionally write one human CDS genomic-position matrix per accepted family
    when `--with-matrices` is used and human is present.
11. Verify the run outputs with `refseq2cds verify`.
12. Optionally align each CDS FASTA with `refseq2cds align` using MAFFT for
    protein MSA and PAL2NAL for codon back-translation.
13. Optionally call target-specific variants with `refseq2cds variants`, using
    the alignment codon maps plus coordinate matrices to write BED rows.
14. Verify alignment and variant output directories with `refseq2cds verify`.

## Install

### Supported Operating Systems

`refseq2cds` is intended to run on macOS and Linux. Windows is not currently
supported. The `refseq2cds download-tools` command detects the operating system
and downloads the matching NCBI Datasets CLI binaries into `./bin`:

```text
macOS -> https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/
Linux -> https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux/
```

The automatic downloader currently supports macOS and Linux only. Even on macOS
or Linux, some environments may require extra setup, especially restricted HPC
systems, containers without `curl`/`unzip`, non-x86_64 CPU architectures, or
networks that block direct access to NCBI FTP/HTTPS endpoints. In those cases,
install NCBI Datasets CLI manually and place `datasets` and `dataformat` in
`./bin` or on your `PATH`.

### Dependencies

Required command-line tools:

| Dependency | Purpose |
|---|---|
| Python `>=3.9` | pipeline runtime |
| `curl` | download NCBI FTP files and NCBI Datasets binaries |
| `gzip` / `unzip` | decompress NCBI inputs |
| NCBI Datasets CLI `datasets` | download RefSeq assembly annotation packages |
| NCBI Datasets CLI `dataformat` | optional NCBI report conversion helper |

Required Python packages are declared in `pyproject.toml` and
`requirements.txt`:

```text
biopython
dendropy
jinja2
networkx
pandas
pyarrow
```

Optional alignment tools for `refseq2cds align --mode mafft-pal2nal`:

| Dependency | Purpose |
|---|---|
| MAFFT | protein multiple sequence alignment |
| PAL2NAL (`pal2nal.pl`) | convert protein MSA + CDS FASTA into codon alignment |
| Perl | runtime for `pal2nal.pl` when installed as a Perl script |

Recommended install with conda/mamba:

```bash
conda install -c bioconda -c conda-forge mafft pal2nal
```

On macOS, MAFFT is also available through Homebrew:

```bash
brew install mafft
```

Homebrew does not always provide PAL2NAL, so use conda/mamba, a system package
manager such as `apt install pal2nal` on Debian/Ubuntu, or pass an explicit
script path with `--pal2nal /path/to/pal2nal.pl`.

Disk requirements depend on the number of input GCF assemblies. The default
14-assembly run currently uses roughly tens of GB because raw NCBI packages and
normalized indexes are kept for reproducibility.

### 1. Clone

```bash
git clone https://github.com/chulbioinfo/refseq2cds.git
cd refseq2cds
```

### 2. Create Python Environment

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

Alternatively, without editable install:

```bash
python -m pip install -r requirements.txt
```

In that mode, run commands as `python refseq2cds.py ...` instead of
`refseq2cds ...`.

Conda/mamba environment:

```bash
mamba env create -f environment.yml
conda activate refseq2cds
```

The environment file installs Python dependencies plus common command-line
tools: NCBI Datasets CLI, MAFFT, and PAL2NAL.

### 3. Download NCBI Datasets CLI

The package can place `datasets` and `dataformat` in `./bin`:

```bash
refseq2cds download-tools
```

or, if running directly from the repository:

```bash
python refseq2cds.py download-tools
```

The workflow was tested with NCBI Datasets CLI `18.25.1`.

## Quick Start

### Offline Smoke Test

Before downloading anything from NCBI, verify the installation with the
packaged mini dataset:

```bash
refseq2cds test --example mini
```

This should finish in seconds. The fixture contains three RefSeq-like species,
one accepted positive-strand family, one accepted negative-strand family, a
missing-species rejection, a paralog-like rejection, and an internal-stop CDS
rejection. It checks graph filtering, CDS QC, FASTA writing, and human
coordinate matrix generation.

For development:

```bash
python -m pip install -e .[test]
pytest -q
```

### Full RefSeq Run

Run the full CDS + matrix workflow with the default 14-species manifest:

```bash
refseq2cds run --steps all --download-tools --with-matrices
```

If you installed with conda/mamba, `ncbi-datasets-cli` is already in the
environment, so `--download-tools` is optional:

```bash
refseq2cds run --steps all --with-matrices
```

By default, gene symbols and FASTA filenames come from human taxid `9606`.
Choose a different reference species with `--reference-taxid`:

```bash
refseq2cds run --steps all --reference-taxid 9598 --with-matrices
```

In that example, singleton family IDs and FASTA filenames use chimpanzee gene
symbols. Human CDS genomic matrices are still generated if human taxid `9606`
is present in the manifest.

Repository-local equivalent:

```bash
./run_cds_pipeline.sh
```

This helper script creates `.venv` when needed, installs Python dependencies,
downloads the matching NCBI Datasets CLI for your OS, and runs the default full
workflow with human matrices.

The first run downloads large NCBI files. Expect tens of GB of local data:

```text
raw/      NCBI bulk files and assembly annotation packages
indexes/  normalized Parquet/TSV indexes
fastas/   final singleton CDS FASTA files
human_cds_matrices/  human CDS-to-genome coordinate matrices
```

Recommended end-to-end command flow:

| Step | Command | Main input | Main output | Used by |
|---|---|---|---|---|
| 1 | `refseq2cds manifest-from-gcf` or `init-manifest` | GCF list or packaged defaults | `config/species_manifest.tsv` | `run`, `verify`, `align` |
| 2 | `refseq2cds run --steps all --with-matrices` | manifest, NCBI Gene bulk files, locked GCF packages | `fastas/`, `human_cds_matrices/`, `reports/` | `verify`, `align`, `variants` |
| 3 | `refseq2cds verify --full --matrix-rows sample` | run output root | JSON pass/fail check | checkpoint before alignment |
| 4 | `refseq2cds align --mode mafft-pal2nal --map-token human` | `fastas/` | `alignments/mafft_pal2nal/` and codon maps | `variants` |
| 5 | `refseq2cds variants` | codon alignments, codon maps, coordinate matrices | `variants/` variant tables and BED files | downstream analysis/UCSC |
| 6 | `refseq2cds verify --alignment-dir ... --variant-dir ...` | completed output directories | JSON pass/fail check | final run audit |

The directory names are intentionally connected: `run` writes `fastas/` and
`human_cds_matrices/`; `align` reads `fastas/` and writes
`alignments/mafft_pal2nal/`; `variants` reads `alignments/mafft_pal2nal/`,
`alignments/mafft_pal2nal/maps/`, and `human_cds_matrices/`.

After FASTA generation, create lightweight tree-free codon alignments:

```bash
refseq2cds align --mode mafft-pal2nal --jobs 4 --threads-per-mafft 2 --map-token human
```

After codon alignment, call target-specific amino acid variants. This example
finds human-specific variants while using human only as the coordinate
reference for BED conversion:

```bash
refseq2cds variants \
  --alignment-dir alignments/mafft_pal2nal \
  --codon-map-dir alignments/mafft_pal2nal/maps \
  --matrix-dir human_cds_matrices \
  --coordinate-reference-token human \
  --target-token human \
  --min-background-informative 5 \
  --force
```

The coordinate reference token is not a privileged biological comparator during
event calling. The command calls codon sites only when target and background
amino acid state sets are mutually exclusive, then maps only events where the
coordinate reference has a real mappable CDS base. v0.1.5 reports coordinateable
substitution-like and gap-involving amino acid differences together as
`bed_event_class=variant`; it does not emit separate insertion/deletion BED
classes. Per-gene BED files are also concatenated and coordinate-sorted into
`variants/merged.bed` for downstream genome-browser and interval workflows. See
`docs/variant_detection.md` for the matrix and BED rules.

### Reference-Gene Query Run

To query one reference gene, for example human `FOXP2`, across the current
manifest:

```bash
refseq2cds run \
  --steps all \
  --orthology-mode reference_gene_1to1_present_species \
  --reference-taxid 9606 \
  --reference-symbol FOXP2 \
  --with-matrices \
  --force
```

If the symbol is ambiguous in the reference species, use the NCBI GeneID
directly:

```bash
refseq2cds run \
  --steps all \
  --orthology-mode reference_gene_1to1_present_species \
  --reference-taxid 9606 \
  --reference-gene-id 93986 \
  --with-matrices \
  --force
```

The default minimum output is two sequences, including the reference sequence.
Change this with `--min-sequences`. Use `--exclude-reference` only if you want
the reference sequence omitted from the FASTA; most comparative analyses should
keep the reference.

With very broad manifests, check disk space before the full assembly-exact run.
Even though reference-gene mode narrows downloads to graph-level 1:1 candidates,
conserved genes may still need many GCF annotation packages. The download-scope
summary is written to `reports/reference_gene_download_scope.json`.

For downstream site mapping, write an alignment-to-CDS codon map for the human
sequence:

```bash
refseq2cds align \
  --mode mafft-pal2nal \
  --jobs 4 \
  --threads-per-mafft 2 \
  --map-token human
```

## Verify Outputs

`refseq2cds verify` is a checkpoint command. It does not regenerate data; it
checks that the files produced by earlier steps are internally consistent.

Fast verification of a run generated with `--with-matrices`:

```bash
refseq2cds verify
```

If outputs are not in the current directory:

```bash
refseq2cds verify --output-root /path/to/run --manifest /path/to/species_manifest.tsv
```

Use `--manifest` with the exact manifest used for the run. FASTA headers are
checked against the manifest species tokens. In strict singleton mode, each
FASTA must contain every manifest token. In reference-gene present-species
mode, the FASTA may contain a validated subset of manifest tokens.

If you ran without `--with-matrices`, skip matrix checks:

```bash
refseq2cds verify --matrix-rows none
```

Full FASTA verification plus sampled matrix row counts:

```bash
refseq2cds verify --full --matrix-rows sample
```

Full matrix row-count verification is slower because it reads all gzipped matrix
files:

```bash
refseq2cds verify --full --matrix-rows full
```

After alignment and variant calling, pass those directories explicitly:

```bash
refseq2cds verify \
  --full \
  --matrix-rows sample \
  --alignment-dir alignments/mafft_pal2nal \
  --variant-dir variants/human_specific
```

What the verify levels mean:

| Option | What is checked | Cost |
|---|---|---|
| default | Python source syntax, required reports, first 200 FASTA files, sampled matrix row counts | fast |
| `--full` | every FASTA file instead of a 200-file sample | moderate for large runs |
| `--matrix-rows none` | skip matrix manifests and matrix row counts | use when no matrices were built |
| `--matrix-rows sample` | count rows for selected matrix files and compare with matrix manifest | fast checkpoint |
| `--matrix-rows full` | count rows in every gzipped matrix file | slow but strongest matrix audit |
| `--alignment-dir PATH` | check alignment `summary.json`, `manifest.tsv`, `failed.tsv`, and codon alignment counts; relative paths are resolved under `--output-root` | fast |
| `--variant-dir PATH` | check variant `summary.json`, `manifest.tsv`, failures, and BED row counts; relative paths are resolved under `--output-root` | fast to moderate |

## Commands

### Initialize Manifest

```bash
refseq2cds init-manifest --force
```

Writes:

```text
config/species_manifest.tsv
```

### Create Manifest From GCF Accessions

Write a file with one RefSeq assembly accession per line:

```text
GCF_009914755.1
GCF_028858775.2
GCF_029289425.2
GCF_029281585.2
GCF_037993035.2-RS_2025_03
GCF_964374335.1/
```

Then generate `config/species_manifest.tsv` automatically and run the pipeline:

```bash
refseq2cds manifest-from-gcf --inputfile gcfs.txt --download-tools --force
refseq2cds run --steps all --force --with-matrices
```

You can also pass accessions directly:

```bash
refseq2cds manifest-from-gcf \
  GCF_009914755.1 GCF_028858775.2 GCF_029289425.2 \
  --force
```

The command queries NCBI Datasets for taxid, organism name, and accession
metadata. Known species from the default manifest receive readable tokens such
as `human`, `chimpanzee`, and `bonobo`; other species receive sanitized
scientific-name tokens. Annotation-release suffixes such as `-RS_2025_03` and
accidental trailing slashes are normalized to the base `GCF_<digits>.<version>`
accession before querying NCBI.
Duplicate taxids are rejected because NCBI Gene orthology is taxid/GeneID based.

Important: changing the manifest changes the biological input set. Use
`refseq2cds run --steps all --force --with-matrices` after generating a new
manifest so old indexes/FASTA files are not reused.

When `--manifest` is omitted, `refseq2cds run` uses
`./config/species_manifest.tsv` if it exists in the current working directory;
otherwise it falls back to the packaged default 14-species manifest.

If you want filenames/symbols from a non-human reference species, pass its
taxid during the run:

```bash
refseq2cds run --steps all --force --reference-taxid 9598 --with-matrices
```

### Run Selected Pipeline Stages

```bash
refseq2cds run --steps preflight,bulk
refseq2cds run --steps assemblies,indexes,edges,singletons
refseq2cds run --steps select,qc,sanity,fastas,report
```

Available stages:

```text
preflight
bulk
assemblies
indexes
edges
singletons
select
qc
sanity
fastas
report
```

### Build Human Coordinate Matrices Only

```bash
refseq2cds build-matrices
```

Overwrite existing matrices:

```bash
refseq2cds build-matrices --force
```

### Build MAFFT + PAL2NAL Codon Alignments

Run all generated FASTA files:

```bash
refseq2cds align --mode mafft-pal2nal --jobs 4 --threads-per-mafft 2
```

Run selected genes:

```bash
refseq2cds align --mode mafft-pal2nal --symbols BRCA1 TP53 A1BG
```

Use an explicit PAL2NAL script path:

```bash
refseq2cds align \
  --mode mafft-pal2nal \
  --pal2nal /path/to/pal2nal.pl
```

Write a human codon map for alignment-site to CDS-site conversion:

```bash
refseq2cds align --mode mafft-pal2nal --map-token human
```

Alignment output:

```text
alignments/mafft_pal2nal/codon/*.codon.fasta
alignments/mafft_pal2nal/aa/*.aa.fasta
alignments/mafft_pal2nal/aa_aligned/*.aa.aln.fasta
alignments/mafft_pal2nal/maps/*.codon_map.tsv.gz
alignments/mafft_pal2nal/manifest.tsv
alignments/mafft_pal2nal/failed.tsv
alignments/mafft_pal2nal/summary.json
```

### Call Target-Specific Variants

`variants` consumes codon alignments, token-specific codon maps, and CDS-to-
genome matrices. It writes per-gene BED tracks plus `merged.bed` for all
coordinateable variants. The default `--alignment-dir` matches the default
`align` output directory:

```bash
refseq2cds variants \
  --alignment-dir alignments/mafft_pal2nal \
  --matrix-dir human_cds_matrices \
  --coordinate-reference-token human \
  --target-token human \
  --min-background-informative 5 \
  --output-dir variants/human_specific \
  --force
```

If `--codon-map-dir` is omitted, the command first looks for
`<alignment-dir>/maps`, which is where `refseq2cds align --map-token human`
writes human codon maps.

### Print Summary

```bash
refseq2cds summary
```

For a separate output directory:

```bash
refseq2cds summary --output-root /path/to/run
```

## Packaging And CI

This repository includes packaging and reviewer-smoke-test scaffolding:

```text
environment.yml              reproducible conda/mamba environment
conda-recipe/meta.yaml       noarch Python conda/Bioconda recipe
.github/workflows/test.yml   GitHub Actions pytest workflow
examples/mini/               offline mini fixture
tests/                       pytest tests for CLI and mini matrix logic
CITATION.cff                 citation metadata
```

The conda recipe is pinned to a GitHub release tarball and SHA256. During
development, update the recipe only after the GitHub tag exists so the final
release archive checksum can be recorded. For example, after tagging v0.1.5,
download the GitHub release tarball, compute its SHA256, update
`conda-recipe/meta.yaml`, and then run the local build checks. To test the
current recipe locally:

```bash
mamba install -c conda-forge -c bioconda conda-build
conda build conda-recipe
```

For a Bioconda submission, copy `conda-recipe/meta.yaml` into
`bioconda-recipes/recipes/refseq2cds/meta.yaml`.

Additional notes:

- `docs/options.md` summarizes commonly used command options.
- `docs/troubleshooting.md` covers common installation and runtime issues.
- `docs/conda_packaging.md` records conda/Bioconda packaging steps.

## Human CDS Genomic Matrix Format

Each matrix is a gzipped TSV. Example:

```text
human_cds_matrices/BRCA1.human_cds_genomic_matrix.tsv.gz
```

Important columns:

| Column | Meaning |
|---|---|
| `family_id` | reference-species symbol + reference-species GeneID family identifier |
| `human_symbol` | original human gene symbol |
| `human_GeneID` | NCBI human GeneID |
| `transcript_accession` | selected human transcript |
| `protein_accession` | selected human protein |
| `cds_nt_pos_1based` | 1-based nucleotide position in normalized human CDS |
| `codon_index_1based` | 1-based codon index in normalized human CDS |
| `codon_offset_1based` | 1, 2, or 3 within codon |
| `cds_base_transcript_orientation` | CDS base in transcript orientation |
| `refseq_accession` | RefSeq genomic sequence accession, e.g. `NC_060941.1` |
| `ucsc_chrom` | UCSC-style chromosome name, e.g. `chr17` |
| `genomic_pos_1based` | 1-based genomic coordinate |
| `bed_start_0based` | BED start coordinate |
| `bed_end_0based` | BED end coordinate |
| `strand` | human CDS strand |
| `gcf_accession` | human assembly accession |
| `annotation_release_id` | RefSeq annotation release ID |

The per-matrix rows describe human CDS nucleotide positions. The accompanying
`human_cds_matrices/manifest.tsv` also records the reference naming fields
(`reference_symbol`, `reference_GeneID`, `reference_taxid`, and
`reference_token`) so non-human reference-species runs remain traceable.

To create BED rows for a codon site, filter by `codon_index_1based` and emit:

```text
ucsc_chrom  bed_start_0based  bed_end_0based  name  score  strand
```

For a nucleotide site, filter by `cds_nt_pos_1based`.

## Example: Inspect BRCA1 Matrix

```bash
gzip -dc human_cds_matrices/BRCA1.human_cds_genomic_matrix.tsv.gz | head
```

## Repository Structure

```text
refseq2cds.py
workflow/scripts/run_cds_pipeline.py
workflow/scripts/build_human_cds_position_matrices.py
workflow/scripts/align_mafft_pal2nal.py
workflow/scripts/call_target_specific_variants.py
config/species_manifest.tsv
examples/mini/
tests/
docs/
conda-recipe/meta.yaml
environment.yml
requirements.txt
pyproject.toml
run_cds_pipeline.sh
```

Large generated directories are intentionally ignored by `.gitignore`:

```text
raw/
indexes/
orthology/
selection/
qc/
fastas/
human_cds_matrices/
reports/
alignments/
```

## Notes And Caveats

- `gene_orthologs.gz` is treated as a pairwise edge list, not a ready-made
  orthogroup table.
- Strict singleton means the induced locked-taxid graph component has exactly one
  GeneID per species and exactly one gene per locked assembly.
- CDS provenance is assembly-exact: selected CDS records must come from the
  frozen `GCF_*` annotation package.
- Human genomic positions use whichever human `GCF_*` assembly is present in
  `config/species_manifest.tsv`. For the default manifest this is
  `GCF_009914755.1` / T2T-CHM13v2.0, not GRCh38.
- Terminal stop codons are removed before FASTA/matrix output and recorded in
  metadata.
- `refseq2cds align --mode mafft-pal2nal` is a lightweight tree-free alignment
  mode. It is useful for fast baseline codon alignments and site-index mapping,
  but it is not a replacement for phylogeny-aware PAGAN2/PRANK analyses when a
  guide tree is central to the study design.
