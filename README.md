# refseq2cds

RefSeq 1:1 orthologs to CDS FASTA files.

Version: `0.1.2`

Author: Chul Lee (chul.bioinfo@gmail.com)

`refseq2cds` is an assembly-exact NCBI RefSeq singleton ortholog CDS FASTA
builder.

This repository builds strict N-way singleton ortholog CDS FASTA files from
NCBI RefSeq assembly annotations and NCBI Gene FTP orthology edges. It also
creates human CDS-position to genomic-position matrices so codon/site-level
results can later be converted to UCSC Genome Browser BED coordinates. It can
also create lightweight tree-free codon alignments with a MAFFT protein MSA
followed by PAL2NAL back-translation.

`refseq2cds` is not a de novo orthology inference program and it is not just a
CDS downloader. It treats NCBI Gene orthology records as input evidence, applies
a strict locked-species singleton filter, selects CDS records from exact RefSeq
`GCF_*` assembly annotation packages, records rejection reasons, and emits
auditable FASTA, metadata, report, and coordinate-matrix outputs.

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

## What It Produces

Main CDS outputs:

```text
fastas/{REFERENCE_SYMBOL}.fasta
fastas/{REFERENCE_SYMBOL}.meta.tsv
fastas/manifest.tsv
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

Each FASTA contains exactly one CDS sequence per locked species. Headers are the locked species
tokens, such as:

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
5. Keep only connected components with exactly one GeneID
   per locked taxid.
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
11. Optionally align each CDS FASTA with `refseq2cds align` using MAFFT for
    protein MSA and PAL2NAL for codon back-translation.

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

After FASTA generation, create lightweight tree-free codon alignments:

```bash
refseq2cds align --mode mafft-pal2nal --jobs 4 --threads-per-mafft 2
```

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

Fast verification:

```bash
refseq2cds verify
```

If outputs are not in the current directory:

```bash
refseq2cds verify --output-root /path/to/run --manifest /path/to/species_manifest.tsv
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

The current generated dataset was verified with:

```text
FASTA files: 13,853
Human CDS genomic matrices: 13,853
Matrix failures: 0
Total matrix rows: 25,337,835
```

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
scientific-name tokens. Duplicate taxids are rejected because NCBI Gene
orthology is taxid/GeneID based.

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

The v0.1.2 conda recipe is pinned to the GitHub release tarball and SHA256. To
test it locally:

```bash
mamba install -c conda-forge -c bioconda conda-build conda-verify
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
