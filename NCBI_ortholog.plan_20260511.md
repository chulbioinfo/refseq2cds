# NCBI Ortholog Assembly-Exact Pipeline Plan

Date: 2026-05-11

## 1. Goal

Build a reproducible pipeline that creates one CDS FASTA file per strict 14-way singleton ortholog family, using:

- NCBI `gene_orthologs.gz` as the only orthology source
- 14 frozen RefSeq assembly accessions (`GCF_*`) as the exact CDS/protein/RNA/GFF provenance
- user-provided ROADIES species tree rooted with `Philippine_flying_lemur` as outgroup
- PAGAN2 as the primary tree-guided codon-aware aligner
- PRANK as fallback, MACSE as diagnostic/QC fallback

The central output is:

```text
fastas/{HUMAN_SYMBOL}.fasta
alignments/final/{HUMAN_SYMBOL}.fas
reports/summary.html
```

Each `fastas/{HUMAN_SYMBOL}.fasta` must contain exactly 14 CDS sequences. FASTA headers must be exactly the locked species tokens, for example `>human`, `>chimpanzee`, `>Philippine_flying_lemur`.

## 2. Locked Policy Decisions

| Item | Decision |
|---|---|
| Orthology source | NCBI Gene FTP `gene_orthologs.gz` only |
| CDS provenance | Assembly-exact mode: download annotation files by frozen `GCF` accession |
| NCBI Datasets group-id | Not used for singleton definition |
| Species set | 14 frozen GCF assemblies |
| Transcript/CDS selection | Assembly-derived CDS only; prefer selected/complete/longest coding isoform |
| Human paralog/many:1 | Reject entire family |
| Lineage expansion | Reject naturally through singleton graph filters |
| MT genes | Exclude before singleton finalization |
| Autosome/sex chromosome | Include autosome and sex chromosome genes |
| Trimming | No alignment trimming; terminal stop removal is CDS normalization, not trimming |
| External validation | NCBI-internal only |
| Refresh | No auto-refresh; all NCBI bulk and assembly packages are timestamped and checksummed |

## 3. Frozen Species Manifest

Create `config/species_manifest.tsv` with these required columns:

```text
token	taxid	scientific_name	common_name	gcf_accession	outgroup
```

| token | taxid | scientific_name | common_name | gcf_accession | outgroup |
|---|---:|---|---|---|---|
| human | 9606 | Homo sapiens | human | GCF_009914755.1 | false |
| chimpanzee | 9598 | Pan troglodytes | chimpanzee | GCF_028858775.2 | false |
| bonobo | 9597 | Pan paniscus | bonobo | GCF_029289425.2 | false |
| gorilla | 9595 | Gorilla gorilla | gorilla | GCF_029281585.2 | false |
| Sumatran_orangutan | 9601 | Pongo abelii | Sumatran orangutan | GCF_028885655.2 | false |
| Bornean_orangutan | 9600 | Pongo pygmaeus | Bornean orangutan | GCF_028885625.2 | false |
| siamang_gibbon | 9590 | Symphalangus syndactylus | siamang gibbon | GCF_028878055.3 | false |
| crab-eating_macaque | 9541 | Macaca fascicularis | crab-eating macaque | GCF_037993035.2 | false |
| pig-tailed_macaque | 9545 | Macaca nemestrina | pig-tailed macaque | GCF_043159975.1 | false |
| common_marmoset | 9483 | Callithrix jacchus | common marmoset | GCF_049354715.1 | false |
| Bolivian_squirrel_monkey | 27679 | Saimiri boliviensis | Bolivian squirrel monkey | GCF_048565385.1 | false |
| sunda_slow_loris | 9470 | Nycticebus coucang | sunda slow loris | GCF_027406575.1 | false |
| ring-tailed_lemur | 9447 | Lemur catta | ring-tailed lemur | GCF_020740605.2 | false |
| Philippine_flying_lemur | 110931 | Cynocephalus volans | Philippine flying lemur | GCF_027409185.1 | true |

Notes:

- `token` is the only allowed FASTA header and tree leaf label.
- `GCF_037993035.2-RS_2025_03` style strings should not be used as assembly accessions. Store the assembly accession as `GCF_037993035.2`; store annotation release IDs separately.

## 4. Required Tools And Packages

### 4.1 System Tools

| Tool | Purpose |
|---|---|
| `curl` | Download NCBI Gene FTP bulk files |
| `sha256sum` or `shasum -a 256` | Checksum all frozen inputs |
| `gzip`, `unzip` | Decompress NCBI files |
| `jq` | Inspect NCBI JSON/JSONL metadata |
| `datasets` | NCBI Datasets CLI |
| `dataformat` | Convert NCBI Datasets reports when useful |
| `snakemake` | Workflow orchestration |
| `apptainer` or `docker` | Run PAGAN2 container |
| `prank` | Alignment fallback |
| `macse` | Frameshift/stop-codon diagnostic fallback |

### 4.2 Python Packages

Use Python 3.11 or newer.

| Package | Purpose |
|---|---|
| `pandas` | Tabular processing |
| `pyarrow` | Parquet output |
| `biopython` | FASTA/GFF sequence handling, translation QC |
| `gffutils` or custom GFF3 parser | Assembly GFF3 indexing |
| `dendropy` | Tree validation |
| `networkx` | Orthology graph singleton checks |
| `pandera` | Optional schema validation |
| `jinja2` | HTML reports |
| `matplotlib` or `plotly` | QC plots |

### 4.3 Conda Environments

`envs/ncbi.yml`

```yaml
name: ncbi-ortholog-ncbi
channels: [conda-forge, bioconda]
dependencies:
  - python>=3.11
  - ncbi-datasets-cli>=16
  - snakemake
  - pandas
  - pyarrow
  - biopython
  - gffutils
  - dendropy
  - networkx
  - pandera
  - jinja2
  - jq
```

`envs/alignment.yml`

```yaml
name: ncbi-ortholog-align
channels: [conda-forge, bioconda]
dependencies:
  - python>=3.11
  - biopython
  - prank
  - macse
```

PAGAN2 should be run via container:

```bash
apptainer build containers/pagan2.sif docker://ariloytynoja/pagan2
```

Wrapper: `bin/pagan2-run`

```bash
#!/usr/bin/env bash
set -euo pipefail
apptainer run -B "$PWD":/data containers/pagan2.sif pagan2 "$@"
```

## 5. Directory Layout

```text
NCBI_ortholog/
├── config/
│   ├── species_manifest.tsv
│   └── pipeline_config.yaml
├── raw/
│   ├── ncbi_bulk/
│   │   ├── gene_orthologs.YYYYMMDD.gz
│   │   ├── gene_info.YYYYMMDD.gz
│   │   └── MANIFEST.checksums
│   └── assembly_packages/
│       └── {token}/
├── indexes/
│   ├── assembly_index.parquet
│   ├── gene_index.parquet
│   ├── transcript_index.parquet
│   ├── cds_index.parquet
│   └── protein_index.parquet
├── orthology/
│   ├── ortholog_edges.parquet
│   ├── candidate_components.parquet
│   ├── strict_singleton.parquet
│   └── rejected.parquet
├── selection/
│   ├── representative_cds.parquet
│   └── selection_audit.parquet
├── qc/
│   ├── cds_qc.parquet
│   ├── family_sanity.parquet
│   └── alignment_qc.parquet
├── fastas/
│   ├── {HUMAN_SYMBOL}.fasta
│   └── {HUMAN_SYMBOL}.meta.tsv
├── trees/
│   ├── species_tree.rooted.binary.nwk
│   └── species_tree.rooted.binary.pagan_safe.nwk
├── alignments/
│   ├── pagan2/
│   ├── prank/
│   └── final/
├── reports/
│   ├── preflight.json
│   ├── provenance.html
│   └── summary.html
├── workflow/
│   ├── Snakefile
│   └── scripts/
├── envs/
├── containers/
└── bin/
```

## 6. Stage-By-Stage Work Plan

### Stage 0. Configuration And Preflight

| Field | Detail |
|---|---|
| Goal | Validate manifest, tools, tree, and GCF availability before large downloads |
| Input | `config/species_manifest.tsv`, `trees/species_tree.rooted.binary.nwk` |
| Tools | Python, `datasets summary genome accession`, `dendropy`, `jq` |
| Packages | `pandas`, `dendropy`, `pandera` |
| Output | `reports/preflight.json`, `reports/gcf_release_ids.tsv`, `reports/chromosome_availability.tsv` |

Checks:

- 14 unique `token` values
- 14 unique `taxid` values
- 14 active `GCF` accessions
- `token` contains no whitespace
- exactly one outgroup, `Philippine_flying_lemur`
- tree leaf set equals manifest token set
- tree is rooted and binary
- first split is `Philippine_flying_lemur` versus the remaining 13 species
- record assembly accession, assembly name, annotation release ID, and annotation release date
- record chromosome/sequence availability, especially MT, X, and Y presence

Important implementation note:

- Parse the manifest as TSV and read only the `token` column for tree validation. Do not use `open(...).read().split()`, because that will mix tokens, taxids, species names, and GCF accessions.

Success metric:

- `reports/preflight.json` has `status: pass`.

### Stage 1. Freeze NCBI Gene Bulk Files

| Field | Detail |
|---|---|
| Goal | Freeze orthology and gene metadata snapshot |
| Input | NCBI Gene FTP |
| Tools | `curl`, `gzip`, `sha256sum` |
| Packages | none |
| Output | `raw/ncbi_bulk/gene_orthologs.YYYYMMDD.gz`, `raw/ncbi_bulk/gene_info.YYYYMMDD.gz`, `raw/ncbi_bulk/MANIFEST.checksums` |

Actions:

```bash
SNAPSHOT_DATE=$(date +%Y%m%d)
mkdir -p raw/ncbi_bulk
curl -L -o raw/ncbi_bulk/gene_orthologs.${SNAPSHOT_DATE}.gz \
  https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz
curl -L -o raw/ncbi_bulk/gene_info.${SNAPSHOT_DATE}.gz \
  https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
(cd raw/ncbi_bulk && sha256sum *.gz > MANIFEST.checksums)
```

Why `gene_info.gz` is needed:

- `gene_orthologs.gz` has only pairwise orthology edges.
- `gene_info.gz` provides `type_of_gene`, `chromosome`, and current symbols for MT/protein-coding filters and human symbol naming.

Success metric:

- Both files decompress successfully.
- Header of `gene_orthologs.gz` is:

```text
#tax_id	GeneID	relationship	Other_tax_id	Other_GeneID
```

### Stage 2. Download Assembly-Exact Annotation Packages

| Field | Detail |
|---|---|
| Goal | Download CDS, protein, RNA, GFF3, genome sequence reports for exactly the frozen `GCF` assemblies |
| Input | `config/species_manifest.tsv` |
| Tools | `datasets download genome accession`, `unzip`, `sha256sum` |
| Packages | NCBI Datasets CLI |
| Output | `raw/assembly_packages/{token}/` |

For each species:

```bash
datasets download genome accession "${GCF}" \
  --include genome,cds,protein,rna,gff3,seq-report \
  --filename "raw/assembly_packages/${TOKEN}.zip" \
  --no-progressbar

unzip -q "raw/assembly_packages/${TOKEN}.zip" \
  -d "raw/assembly_packages/${TOKEN}"
```

Implementation requirements:

- Do not hardcode exact filenames such as `cds.fna`; inspect `dataset_catalog.json` and locate files by type.
- Record package checksum.
- Record package `GCF`, annotation release ID, release date, and file list.
- Fail if downloaded assembly accession does not equal manifest `gcf_accession`.

Expected files may include:

- genomic FASTA
- CDS FASTA
- protein FASTA
- RNA FASTA
- GFF3
- sequence report
- assembly data report
- dataset catalog

Success metric:

- 14/14 assembly packages downloaded.
- Each package contains GFF3, CDS FASTA, protein FASTA, RNA FASTA, and metadata reports.

### Stage 3. Build Assembly-Exact Indexes

| Field | Detail |
|---|---|
| Goal | Convert assembly package contents into normalized tables keyed by taxid, GeneID, transcript accession, CDS accession, and protein accession |
| Input | `raw/assembly_packages/{token}/`, `raw/ncbi_bulk/gene_info.*.gz` |
| Tools | Python |
| Packages | `pandas`, `pyarrow`, `biopython`, `gffutils` or a strict GFF3 parser |
| Output | `indexes/assembly_index.parquet`, `indexes/gene_index.parquet`, `indexes/transcript_index.parquet`, `indexes/cds_index.parquet`, `indexes/protein_index.parquet` |

Index schemas:

`assembly_index.parquet`

```text
token, taxid, gcf_accession, assembly_name, annotation_release_id,
annotation_release_date, package_path, gff3_path, cds_fasta_path,
rna_fasta_path, protein_fasta_path, seq_report_path
```

`gene_index.parquet`

```text
token, taxid, GeneID, symbol, chromosome, type_of_gene,
feature_id, gff3_seqid, start, end, strand, source
```

`transcript_index.parquet`

```text
token, taxid, GeneID, transcript_accession, transcript_feature_id,
transcript_type, select_category, transcript_length,
gff3_seqid, start, end, strand, source
```

`cds_index.parquet`

```text
token, taxid, GeneID, transcript_accession, cds_accession,
protein_accession, cds_length, cds_fasta_header, cds_sequence,
is_partial, source
```

`protein_index.parquet`

```text
token, taxid, GeneID, transcript_accession, protein_accession,
protein_length, protein_fasta_header, source
```

Parsing rules:

- Derive GeneID from GFF3 `Dbxref=GeneID:<id>` when available.
- Derive transcript accession from RNA/mRNA feature attributes and CDS parent relationships.
- Derive protein accession from GFF3 `protein_id`, CDS FASTA headers, or protein FASTA headers.
- Prefer GFF3 feature relationships over FASTA header parsing.
- Use FASTA headers only as secondary reconciliation.
- Keep `source = assembly_exact`.

Success metrics:

- All `cds_index` rows map to one taxid and one frozen assembly package.
- No orphan CDS without GeneID unless explicitly logged and excluded.
- No duplicate primary key for `(taxid, GeneID, transcript_accession, cds_accession)`.
- All retained CDS rows have a sequence.

### Stage 4. Build Orthology Edge Graph From `gene_orthologs.gz`

| Field | Detail |
|---|---|
| Goal | Parse pairwise NCBI orthology edges among the 14 taxids |
| Input | `raw/ncbi_bulk/gene_orthologs.*.gz`, `config/species_manifest.tsv` |
| Tools | Python |
| Packages | `pandas`, `pyarrow`, `networkx` |
| Output | `orthology/ortholog_edges.parquet` |

Actions:

- Keep only rows with `relationship == Ortholog`.
- Keep only edges where both `tax_id` and `Other_tax_id` are in the 14 locked taxids.
- Normalize edge direction so `(taxid_a, GeneID_a, taxid_b, GeneID_b)` is deterministic.
- Remove duplicate edges.
- Build an undirected graph over nodes `(taxid, GeneID)`.

Output schema:

```text
taxid_a, GeneID_a, taxid_b, GeneID_b, relationship, source_snapshot
```

Success metric:

- All edges are internal to the 14-species set.
- `relationship` is always `Ortholog`.

### Stage 5. Strict 14-Way Singleton Graph Filter

| Field | Detail |
|---|---|
| Goal | Define strict singleton ortholog families as graph components with exactly one GeneID per species |
| Input | `orthology/ortholog_edges.parquet`, `indexes/gene_index.parquet`, `raw/ncbi_bulk/gene_info.*.gz` |
| Tools | Python |
| Packages | `pandas`, `networkx`, `pyarrow` |
| Output | `orthology/candidate_components.parquet`, `orthology/strict_singleton.parquet`, `orthology/rejected.parquet` |

Filtering logic:

1. Build connected components from the 14-species orthology graph.
2. Keep only components that contain `human` taxid 9606.
3. Reject component if it contains any gene with `chromosome == MT`.
4. Reject component if any gene is not protein-coding according to `gene_info.gz` or assembly GFF3 metadata.
5. Reject component if human GeneID count is not exactly 1.
6. Reject component if any taxid has 0 genes.
7. Reject component if any taxid has more than 1 GeneID.
8. Reject component if total gene count is not exactly 14.
9. Reject component if any component gene is absent from `gene_index.parquet`.
10. Keep component as strict singleton if all checks pass.

Reject reason categories:

```text
missing_human
multiple_human_genes
missing_taxon
multi_gene_per_taxon
component_size_not_14
mt_gene
non_protein_coding
missing_assembly_gene
no_ortholog_edge_within_locked_taxa
```

Output `strict_singleton.parquet` schema:

```text
family_id, human_GeneID, human_symbol,
taxid, token, GeneID, species_symbol, component_id
```

Success metrics:

- Every accepted `family_id` has exactly 14 rows.
- Every accepted `family_id` has exactly 14 unique taxids.
- Every accepted `family_id` has exactly one human GeneID.
- Rejected components have exactly one primary reject reason.

### Stage 6. Select Representative Assembly-Exact CDS

| Field | Detail |
|---|---|
| Goal | Select one CDS per `(family_id, taxid)` from the frozen assembly annotation |
| Input | `orthology/strict_singleton.parquet`, `indexes/transcript_index.parquet`, `indexes/cds_index.parquet`, `indexes/protein_index.parquet` |
| Tools | Python |
| Packages | `pandas`, `pyarrow` |
| Output | `selection/representative_cds.parquet`, `selection/selection_audit.parquet` |

Selection rule:

1. Candidate must come from the frozen assembly package.
2. Candidate must have a CDS sequence.
3. Candidate must have a protein accession or valid coding translation unless documented otherwise.
4. Exclude partial CDS candidates when partial status can be determined.
5. Prefer `select_category` values such as `MANE Select` or `RefSeq Select` when present.
6. Prefer complete protein-coding transcript.
7. Prefer longest protein length.
8. Prefer longest CDS length.
9. Prefer curated accession class (`NM_`) over predicted accession class (`XM_`) only after assembly-exact and coding-completeness criteria.
10. Use lexicographically smallest transcript accession as deterministic final tie-break.

Important distinction:

- Do not select by longest RNA/transcript length alone.
- The intended representative is the longest defensible coding isoform, not the longest UTR-containing transcript.

Output schema:

```text
family_id, human_symbol, taxid, token, GeneID, species_symbol,
transcript_accession, cds_accession, protein_accession,
cds_length, protein_length, select_category, accession_class,
selection_rule_id, source
```

Success metrics:

- Every accepted `family_id` has exactly 14 selected CDS records.
- Every `(family_id, taxid)` has exactly one selected CDS.
- All selected CDS have `source == assembly_exact`.

### Stage 7. CDS Sequence QC And Normalization

| Field | Detail |
|---|---|
| Goal | Ensure all selected CDS sequences are safe for PAGAN2 codon mode |
| Input | `selection/representative_cds.parquet`, `indexes/cds_index.parquet` |
| Tools | Python |
| Packages | `biopython`, `pandas`, `pyarrow` |
| Output | `qc/cds_qc.parquet`, normalized sequences in memory or `selection/representative_cds.normalized.parquet` |

Checks:

- length >= 60 nt
- length divisible by 3 after terminal stop normalization
- valid DNA alphabet; allow only `ACGTN` by default
- `N` fraction <= 1% by default
- no internal stop codon under standard genetic code
- terminal stop codon may be removed and recorded
- optional start codon check: record but do not hard-fail by default

Terminal stop policy:

- If final codon is a stop codon, remove it before alignment.
- Record `had_terminal_stop = true`.
- Record `terminal_stop_removed = true`.
- This is CDS normalization, not alignment trimming.

Family-level policy:

- If any one species fails CDS QC, reject the entire family from FASTA/alignment generation.

Success metrics:

- All families moving forward have 14/14 QC-passed normalized CDS.
- No sequence entering PAGAN2 contains internal stops.

### Stage 8. Family Sanity Check

| Field | Detail |
|---|---|
| Goal | Flag possible residual paralogy, bad annotation, or severe isoform mismatch before alignment |
| Input | QC-passed CDS/protein records |
| Tools | Python; optional fast protein aligner for pilot only |
| Packages | `pandas`, `biopython` |
| Output | `qc/family_sanity.parquet` |

Recommended metrics:

- protein length min, median, max
- protein length max/min ratio
- CDS length max/min ratio
- per-species length z-score
- optional rough protein identity sketch

Initial policy:

- Treat sanity thresholds as `flag`, not hard fail, during MVP.
- After pilot distribution is observed, promote extreme thresholds to hard fail if needed.

Default flags:

```text
protein_length_ratio_gt_1_5
protein_length_ratio_gt_2_0
cds_length_ratio_gt_2_0
species_length_outlier
```

Success metric:

- Every retained family has a sanity status of `pass` or `flag`.

### Stage 9. Write Human-Symbol FASTA Files

| Field | Detail |
|---|---|
| Goal | Create one unaligned CDS FASTA per strict singleton family |
| Input | `selection/representative_cds.normalized.parquet`, `qc/cds_qc.parquet`, `qc/family_sanity.parquet`, `orthology/strict_singleton.parquet` |
| Tools | Python |
| Packages | `biopython`, `pandas` |
| Output | `fastas/{HUMAN_SYMBOL}.fasta`, `fastas/{HUMAN_SYMBOL}.meta.tsv`, `fastas/manifest.tsv` |

Rules:

- Filename is sanitized human symbol.
- Original human symbol is preserved in metadata.
- If sanitized filename collides, append human GeneID: `{SYMBOL}__{human_GeneID}.fasta`.
- FASTA header is exactly `token`, with no whitespace.
- Sequence order follows `species_manifest.tsv`.
- Use normalized CDS sequence with terminal stop removed where applicable.

FASTA example:

```text
>human
ATG...
>chimpanzee
ATG...
>bonobo
ATG...
```

Metadata fields:

```text
family_id, human_symbol, fasta_path, token, taxid, GeneID,
species_symbol, transcript_accession, cds_accession,
protein_accession, cds_length_original, cds_length_normalized,
protein_length, select_category, accession_class,
had_terminal_stop, terminal_stop_removed, gcf_accession,
annotation_release_id
```

Success metrics:

- Every FASTA has exactly 14 records.
- Header set equals manifest token set.
- No filename collisions unresolved.

### Stage 10. Validate User-Provided ROADIES Tree

| Field | Detail |
|---|---|
| Goal | Confirm the tree is acceptable for PAGAN2 and PRANK |
| Input | `trees/species_tree.rooted.binary.nwk`, `config/species_manifest.tsv` |
| Tools | Python |
| Packages | `dendropy` |
| Output | `trees/species_tree.rooted.binary.pagan_safe.nwk`, `reports/tree_validation.json` |

Checks:

- leaf set equals manifest token set
- rooted
- binary
- `Philippine_flying_lemur` is outgroup
- all branch lengths exist
- branch lengths are positive

Normalization:

- If branch length is missing or zero, write a PAGAN-safe copy with epsilon branch length, for example `1e-6`.
- Do not silently reroot or resolve polytomies in the main pipeline. Fail with a clear message if the user-provided tree is not rooted and binary.

Success metric:

- `species_tree.rooted.binary.pagan_safe.nwk` exists and passes validation.

### Stage 11. Tree-Guided Codon-Aware Alignment

| Field | Detail |
|---|---|
| Goal | Align each singleton CDS FASTA using the validated ROADIES species tree |
| Input | `fastas/{HUMAN_SYMBOL}.fasta`, `trees/species_tree.rooted.binary.pagan_safe.nwk` |
| Tools | PAGAN2 container, PRANK, MACSE |
| Packages | `biopython` for wrapper QC |
| Output | `alignments/pagan2/{HUMAN_SYMBOL}.fas`, `alignments/prank/{HUMAN_SYMBOL}.best.fas`, `alignments/final/{HUMAN_SYMBOL}.fas`, `alignments/aligner_log.parquet` |

Recommended implementation:

- Do not model PAGAN2 and PRANK fallback as fragile separate Snakemake graph branches.
- Implement one wrapper script: `workflow/scripts/11_align_one.py`.
- For each family, the wrapper should:
  1. run PAGAN2
  2. run immediate alignment QC
  3. accept PAGAN2 if QC passes
  4. otherwise run PRANK
  5. run immediate alignment QC
  6. accept PRANK if QC passes
  7. otherwise optionally run MACSE diagnostic and mark failed

PAGAN2 command pattern:

```bash
bin/pagan2-run \
  --seqfile /data/fastas/${SYMBOL}.fasta \
  --treefile /data/trees/species_tree.rooted.binary.pagan_safe.nwk \
  --codons \
  --outfile /data/alignments/pagan2/${SYMBOL} \
  --outformat fasta \
  --config-log-file /data/alignments/pagan2/${SYMBOL}.cfg
```

PRANK fallback command pattern:

```bash
prank \
  -d=fastas/${SYMBOL}.fasta \
  -t=trees/species_tree.rooted.binary.pagan_safe.nwk \
  -codon \
  -prunedata \
  -prunetree \
  -o=alignments/prank/${SYMBOL}
```

Important:

- PAGAN/PAGAN2 codon mode assumes first reading frame and does not correct frameshifts.
- Stage 7 QC is therefore mandatory.

Success metrics:

- PAGAN2 success rate reported.
- PRANK fallback success rate reported.
- Every accepted final alignment has exactly 14 sequences.

### Stage 12. Alignment QC

| Field | Detail |
|---|---|
| Goal | Verify final alignments are usable for downstream codon analyses |
| Input | `alignments/final/*.fas`, `trees/species_tree.rooted.binary.pagan_safe.nwk` |
| Tools | Python |
| Packages | `biopython`, `pandas`, `pyarrow` |
| Output | `qc/alignment_qc.parquet` |

Checks:

- sequence count == 14
- header set equals tree leaf set
- alignment length divisible by 3
- no internal stop codons after removing gaps by codon
- per-sequence gap fraction
- per-family gap fraction
- species-specific outlier gaps
- no sequence dropped by aligner

Initial policy:

- gap fraction > 50%: flag
- gap fraction > 80%: fail
- internal stop: fail
- length not divisible by 3: fail

Success metrics:

- Every accepted alignment has `status == pass` or `status == flag`.
- Every failed alignment has a categorized fail reason.

### Stage 13. Final Report And Provenance

| Field | Detail |
|---|---|
| Goal | Summarize all inputs, filters, failures, and final outputs |
| Input | All stage outputs |
| Tools | Python |
| Packages | `pandas`, `jinja2`, `matplotlib` or `plotly` |
| Output | `reports/summary.html`, `reports/provenance.html`, `reports/final_manifest.tsv` |

Report contents:

- run date
- git commit if applicable
- NCBI Datasets CLI version
- PAGAN2 container digest or image ID
- PRANK and MACSE versions
- NCBI Gene FTP snapshot date and checksums
- 14 GCF accessions and annotation release IDs
- number of graph components
- number of strict singleton families
- rejection counts by stage and reason
- selected transcript/CDS accession class distribution
- per-species selected `NM_` versus `XM_` counts
- CDS QC fail counts
- family sanity flags
- PAGAN2/PRANK/fail alignment counts
- final accepted alignment count

Success metric:

- `reports/final_manifest.tsv` links every final alignment back to:
  - human symbol
  - human GeneID
  - family ID
  - species token
  - taxid
  - species GeneID
  - transcript accession
  - CDS accession
  - protein accession
  - GCF accession
  - annotation release ID

## 7. Minimal Snakemake Shape

Use a checkpoint or generated manifest after FASTA creation. Avoid using `glob("fastas/*.fasta")` at DAG parse time before FASTA files exist.

Conceptual flow:

```text
preflight
bulk_freeze
download_assembly_packages
build_indexes
build_orthology_graph
strict_singleton
representative_cds
cds_qc
family_sanity
write_fastas
validate_tree
align_all_from_fasta_manifest
alignment_qc
summary
```

Recommended alignment rule:

```python
rule align_one:
    input:
        fasta="fastas/{symbol}.fasta",
        tree="trees/species_tree.rooted.binary.pagan_safe.nwk"
    output:
        aln="alignments/final/{symbol}.fas",
        status="alignments/final/{symbol}.status.json"
    log:
        "logs/align/{symbol}.log"
    conda:
        "envs/alignment.yml"
    shell:
        """
        python workflow/scripts/11_align_one.py \
          --fasta {input.fasta} \
          --tree {input.tree} \
          --symbol {wildcards.symbol} \
          --out {output.aln} \
          --status {output.status} \
          > {log} 2>&1
        """
```

## 8. Main Risks And Mitigations

| Risk | Why It Matters | Mitigation |
|---|---|---|
| `gene_orthologs.gz` is pairwise, not a final 14-way orthogroup table | Human-anchor aggregation can miss complex graph structure | Build full 14-taxid graph and require component size 14 with one GeneID per taxid |
| Gene metadata and assembly annotation disagree | Gene FTP is current, while GCF annotation is frozen | Use assembly GFF/CDS as sequence truth; use `gene_info.gz` only for supplementary metadata and flag disagreements |
| GFF3 attribute parsing varies by assembly | CDS-transcript-GeneID links may be missed | Prefer formal GFF3 parent relationships; log orphan features; use FASTA headers only as reconciliation |
| Longest transcript may select UTR-heavy RNA | Bad codon alignment target | Select longest coding isoform by protein/CDS length, not RNA length |
| `NM_` priority may conflict with assembly-exact provenance | Curated transcript may not be present in frozen GCF annotation | Apply `NM_ > XM_` only after assembly-exact candidate filtering |
| PAGAN2 does not fix frames | Bad CDS enters alignment and creates invalid codon MSA | Strict Stage 7 CDS QC before alignment |
| Tree labels mismatch FASTA headers | Aligner failure or wrong mapping | Manifest token is single source of truth; Stage 10 hard-fails mismatch |
| Snakemake fallback logic becomes brittle | Failed PAGAN2 jobs can stop workflow prematurely | Use one `align_one.py` wrapper that handles PAGAN2, QC, PRANK fallback, and status logging |

## 9. MVP Milestones

| Milestone | Scope | Exit Criteria |
|---|---|---|
| MVP 0 | Manifest, preflight, NCBI bulk freeze | `reports/preflight.json` pass; checksums recorded |
| MVP 1 | Assembly package download and indexes | 14/14 GCF packages indexed; CDS/GFF/protein/RNA tables built |
| MVP 2 | Orthology graph and strict singleton filter | `orthology/strict_singleton.parquet` and rejection summary produced |
| MVP 3 | Representative CDS and FASTA generation | `fastas/manifest.tsv` created; every FASTA has 14 sequences |
| V1 | PAGAN2/PRANK alignment and QC | final alignments and `qc/alignment_qc.parquet` created |
| V1 report | Full provenance | `reports/summary.html` and `reports/final_manifest.tsv` complete |

## 10. Definition Of Done

The pipeline is considered complete when a clean run from frozen inputs produces:

```text
orthology/strict_singleton.parquet
selection/representative_cds.parquet
qc/cds_qc.parquet
fastas/manifest.tsv
alignments/final/*.fas
qc/alignment_qc.parquet
reports/summary.html
reports/provenance.html
reports/final_manifest.tsv
```

All final alignments must satisfy:

- exactly 14 sequences
- headers equal the 14 locked species tokens
- alignment length divisible by 3
- no internal stop codons
- full provenance back to `gene_orthologs.gz`, `GeneID`, `GCF`, transcript accession, CDS accession, and protein accession

