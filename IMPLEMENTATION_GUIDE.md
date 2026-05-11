# 실무 지침서: NCBI gene_orthologs.gz 기반 14종 strict singleton CDS + codon alignment

**작성일**: 2026-05-10
**목적**: 14 영장류(13 primates + Philippine flying lemur outgroup)에 대해 NCBI 단독 의존(strict, fast) 1:1 singleton CDS 추출 → ROADIES species tree guide → PAGAN2 codon-aware alignment 자동화

---

## 목차

- [0. 프로젝트 개요](#0-프로젝트-개요)
- [1. 디렉토리 구조](#1-디렉토리-구조)
- [2. 환경 설정](#2-환경-설정)
- [3. 단계별 실행 계획](#3-단계별-실행-계획)
  - [Stage 0: Pre-flight 검증](#stage-0-pre-flight-검증)
  - [Stage 1: NCBI bulk 데이터 freeze](#stage-1-ncbi-bulk-데이터-freeze)
  - [Stage 2: 종별 gene package 다운로드](#stage-2-종별-gene-package-다운로드)
  - [Stage 3: 데이터 정규화](#stage-3-데이터-정규화-parquet-index)
  - [Stage 4: Orthology backbone](#stage-4-orthology-backbone-from-gene_orthologsgz)
  - [Stage 5: Strict singleton + human paralog 배제 + MT 배제](#stage-5-strict-singleton--human-paralog-배제--mt-배제)
  - [Stage 6: 대표 CDS 선택](#stage-6-종별-대표-cds-선택-nm_-우선-xm_-fallback)
  - [Stage 7: CDS QC](#stage-7-cds-무결성-qc)
  - [Stage 8: 종간 sanity check](#stage-8-종간-sanity-check)
  - [Stage 9: FASTA 생성](#stage-9-symbolfasta-생성)
  - [Stage 10: Tree 검증](#stage-10-roadies-tree-검증-사용자-제공)
  - [Stage 11: Codon alignment](#stage-11-pagan2-codon-aware-alignment)
  - [Stage 12: Alignment QC](#stage-12-alignment-qc)
  - [Stage 13: 최종 리포트](#stage-13-최종-리포트--provenance)
- [4. Snakemake workflow 스켈레톤](#4-snakemake-workflow-스켈레톤)
- [5. 리스크 매트릭스](#5-리스크-매트릭스)
- [6. 즉시 시작할 수 있는 명령어](#6-즉시-시작할-수-있는-명령어-3시간-내-mvp-0-도달)
- [7. 다음 작업 권장 순서](#7-다음-작업-권장-순서)

---

## 0. 프로젝트 개요

### 0.1 Lock된 정책 결정 (9개)

| # | 항목 | 결정 |
|---|---|---|
| 1 | Orthology 소스 | `gene_orthologs.gz` 단독 (NCBI Datasets group-id 사용 안 함) |
| 2 | Transcript | NM_ 우선, 없으면 XM_ |
| 3 | Tree | 사용자가 ROADIES + flying lemur outgroup으로 제공 |
| 4 | MT | 사전 배제 |
| 5 | Human paralog | `gene_orthologs.gz`에서 many:1 발생하는 모든 family 배제 |
| 6 | Lineage expansion | singleton 필터로 자연 reject. autosome + sex chr 포함 |
| 7 | Trimming | 없음 (untrimmed alignment만) |
| 8 | External validation | NCBI 내부만 (OrthoMaM/OMA/Ensembl 사용 안 함) |
| 9 | Refresh | 14 GCF로 freeze, no auto-refresh |

### 0.2 14종 freeze (file name token 포함)

| token (header/leaf) | taxid | GCF | 비고 |
|---|---|---|---|
| human | 9606 | GCF_009914755.1 | T2T-CHM13v2.0 (비표준 reference) |
| chimpanzee | 9598 | GCF_028858775.2 | |
| bonobo | 9597 | GCF_029289425.2 | |
| gorilla | 9595 | GCF_029281585.2 | |
| Sumatran_orangutan | 9601 | GCF_028885655.2 | |
| Bornean_orangutan | 9600 | GCF_028885625.2 | |
| siamang_gibbon | 9590 | GCF_028878055.3 | |
| crab-eating_macaque | 9541 | GCF_037993035.2 | RS_YYYY_MM suffix 주의 |
| pig-tailed_macaque | 9545 | GCF_043159975.1 | |
| common_marmoset | 9483 | GCF_049354715.1 | |
| Bolivian_squirrel_monkey | 27679 | GCF_048565385.1 | |
| sunda_slow_loris | 9470 | GCF_027406575.1 | |
| ring-tailed_lemur | 9447 | GCF_020740605.2 | |
| Philippine_flying_lemur | 110931 | GCF_027409185.1 | **outgroup** |

---

## 1. 디렉토리 구조

```
NCBI_ortholog/
├── config/
│   ├── species_manifest.tsv          # 14종 metadata (locked)
│   └── pipeline_config.yaml          # threshold/path/version
├── raw/
│   ├── ncbi_bulk/
│   │   ├── gene_orthologs.<date>.gz  # frozen snapshot
│   │   ├── gene_info.<date>.gz       # MT/chromosome/gene_type 필터용
│   │   └── MANIFEST.checksums
│   ├── gene_packages/
│   │   └── {taxid}/                  # datasets download gene taxon 산출물
│   └── genomes/                      # ROADIES용 (사용자가 별도 처리)
├── indexes/
│   ├── gene_index.parquet            # 14종 통합 gene metadata
│   ├── product_index.parquet         # transcript/CDS/protein metadata
│   └── cds_index.parquet
├── orthology/
│   ├── ortholog_candidates.parquet
│   ├── strict_singleton.parquet
│   └── rejected.parquet              # 사유 분류된 reject log
├── selection/
│   ├── representative_cds.parquet
│   └── selection_audit.parquet
├── qc/
│   ├── cds_qc.parquet
│   ├── family_sanity.parquet
│   └── alignment_qc.parquet
├── fastas/
│   ├── {HUMAN_SYMBOL}.fasta          # 14 sequence per file
│   └── {HUMAN_SYMBOL}.meta.tsv
├── trees/
│   └── species_tree.rooted.binary.nwk  # 사용자 제공 (검증 후 보관)
├── alignments/
│   └── pagan2/
│       ├── {HUMAN_SYMBOL}.fas
│       └── {HUMAN_SYMBOL}.log
├── reports/
│   ├── preflight.json
│   ├── provenance.html
│   └── summary.html
├── workflow/
│   ├── Snakefile
│   └── scripts/
├── containers/
│   └── pagan2.sif                    # Apptainer image
└── bin/
    └── pagan2-run                    # wrapper
```

---

## 2. 환경 설정

### 2.1 Conda 환경 (2개로 분리)

```yaml
# envs/ncbi.yml
name: ncbi
channels: [conda-forge, bioconda]
dependencies:
  - ncbi-datasets-cli>=16
  - python>=3.11
  - pandas
  - pyarrow
  - biopython
  - dendropy           # tree 검증용
  - requests
  - jq
```

```yaml
# envs/alignment.yml
name: align
channels: [conda-forge, bioconda]
dependencies:
  - prank             # PAGAN2 fallback
  - macse>=2          # frame 분기용 (필요 시)
  - python>=3.11
  - biopython
```

### 2.2 PAGAN2 컨테이너 (Docker 또는 Apptainer)

```bash
# Docker 환경
docker pull ariloytynoja/pagan2
docker tag ariloytynoja/pagan2 pagan2

# HPC (Apptainer/Singularity) 환경
apptainer build containers/pagan2.sif docker://ariloytynoja/pagan2
```

```bash
# bin/pagan2-run (Apptainer 가정, Docker는 docker run으로 치환)
#!/usr/bin/env bash
set -euo pipefail
apptainer run -B "$PWD":/data containers/pagan2.sif pagan2 "$@"
```

---

## 3. 단계별 실행 계획

### Stage 0: Pre-flight 검증

| 항목 | 내용 |
|---|---|
| **목표** | 본격 다운로드 전 14 GCF 메타데이터 정합성, NCBI API 가용성, 사용자 제공 tree leaf 일치성 검증 |
| **인풋** | `config/species_manifest.tsv`, 사용자 제공 `trees/species_tree.rooted.binary.nwk` |
| **도구** | `datasets summary genome accession`, dendropy, Python |
| **액션** | (a) 각 GCF의 `assembly_info.annotation_info.release_id` 추출 → `gcf_release_ids.tsv`<br>(b) 각 GCF의 chromosome 목록·`has_Y` 추출 → `chromosome_availability.tsv`<br>(c) tree leaf set이 manifest의 14 token과 정확히 일치하는지 확인<br>(d) tree가 rooted binary인지, flying lemur가 outgroup 위치인지 검증 |
| **아웃풋** | `reports/preflight.json` (모든 검사 결과), 실패 시 즉시 중단 |
| **검증** | 14/14 GCF active, 14/14 release_id 기록, tree leaf == species token, rooted=true, binary=true |

```bash
# (a) annotation release ID 수집
while IFS=$'\t' read -r token taxid gcf; do
  datasets summary genome accession "$gcf" --as-json-lines \
    | jq -r '[.accession, .annotation_info.release_id // "NA", .annotation_info.release_date // "NA"] | @tsv'
done < <(tail -n +2 config/species_manifest.tsv | cut -f1,2,3) \
  > reports/gcf_release_ids.tsv
```

```python
# (c)(d) tree 검증 (dendropy)
import dendropy
tree = dendropy.Tree.get(path="trees/species_tree.rooted.binary.nwk", schema="newick")
manifest_tokens = {line.split("\t")[0] for line in
                   open("config/species_manifest.tsv").readlines()[1:]}
leaf_tokens = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
assert leaf_tokens == manifest_tokens, f"mismatch: {leaf_tokens ^ manifest_tokens}"
assert tree.is_rooted
assert all(len(node.child_nodes()) in (0, 2) for node in tree.preorder_node_iter())
```

---

### Stage 1: NCBI bulk 데이터 freeze

| 항목 | 내용 |
|---|---|
| **목표** | `gene_orthologs.gz`와 `gene_info.gz`(MT/chromosome/gene_type 필터용)를 단일 timestamp로 freeze |
| **인풋** | 없음 (외부 NCBI FTP) |
| **도구** | `curl`, `sha256sum` |
| **액션** | NCBI Gene FTP에서 `gene_orthologs.gz`, `gene_info.gz` 다운로드 후 SHA256 기록 |
| **아웃풋** | `raw/ncbi_bulk/gene_orthologs.<YYYYMMDD>.gz`, `raw/ncbi_bulk/gene_info.<YYYYMMDD>.gz`, `MANIFEST.checksums` |
| **검증** | 파일 크기 정상, SHA256 기록, 다운로드 timestamp 기록 |

```bash
SNAPSHOT_DATE=$(date +%Y%m%d)
mkdir -p raw/ncbi_bulk
for f in gene_orthologs.gz gene_info.gz; do
  curl -sS -o "raw/ncbi_bulk/${f%.gz}.${SNAPSHOT_DATE}.gz" \
    "https://ftp.ncbi.nlm.nih.gov/gene/DATA/$f"
done
(cd raw/ncbi_bulk && sha256sum *.gz > MANIFEST.checksums)
```

**`gene_info.gz`가 필요한 이유**: MT 배제(정책 #4)와 `type_of_gene == 'protein-coding'` 필터를 위해 chromosome·gene_type 컬럼이 필요한데, `gene_orthologs.gz`에는 이 정보가 없습니다.

---

### Stage 2: 종별 gene package 다운로드

| 항목 | 내용 |
|---|---|
| **목표** | 14종 각각의 gene/CDS/product report 패키지 수집 |
| **인풋** | `species_manifest.tsv` (taxid 14개) |
| **도구** | NCBI Datasets CLI (`datasets download gene taxon ... --include cds,product-report`) |
| **액션** | taxid별로 gene package 다운로드, 압축 해제, 재시도 backoff 적용 |
| **아웃풋** | `raw/gene_packages/{taxid}/ncbi_dataset/data/cds.fna`, `product_report.jsonl`, `gene_report.jsonl` 등 |
| **검증** | 14/14 taxid 성공, 각 package에 cds.fna 존재, product_report 비어있지 않음 |

```bash
mkdir -p raw/gene_packages
while IFS=$'\t' read -r token taxid gcf; do
  out="raw/gene_packages/${taxid}"
  [ -d "$out" ] && continue
  for attempt in 1 2 3; do
    datasets download gene taxon "$taxid" \
      --include cds,product-report,gene \
      --filename "raw/gene_packages/${taxid}.zip" \
      --no-progressbar && break
    sleep $((10 * attempt))
  done
  unzip -q "raw/gene_packages/${taxid}.zip" -d "$out"
  rm "raw/gene_packages/${taxid}.zip"
done < <(tail -n +2 config/species_manifest.tsv | cut -f1,2,3)
```

**주의**: GCF↔annotation release 정합성에서 release_id가 manifest와 다르면 fail시키거나, GCF별 annotation 직접 다운로드 (genome accession 기반) 대안을 사용.

---

### Stage 3: 데이터 정규화 (parquet index)

| 항목 | 내용 |
|---|---|
| **목표** | 14종 raw 데이터를 join 가능한 정형 테이블로 변환 |
| **인풋** | `raw/gene_packages/{taxid}/*`, `raw/ncbi_bulk/gene_info.*.gz` |
| **도구** | Python (pandas + pyarrow), `dataformat tsv gene`, `dataformat tsv gene-product` |
| **액션** | (a) `gene_index.parquet`: 14종 통합 `(taxid, GeneID, symbol, chromosome, type_of_gene, gene_biotype)` — `gene_info.gz`와 gene_report 조인<br>(b) `product_index.parquet`: `(GeneID, transcript_accession, transcript_length, cds_accession, protein_accession, protein_length, transcript_select_category, accession_class)` — `accession_class`는 `NM_`/`XM_`/`NR_`/`XR_` 접두사로 도출<br>(c) `cds_index.parquet`: `(cds_accession, taxid, GeneID, sequence_length, file_path)` |
| **아웃풋** | `indexes/{gene,product,cds}_index.parquet` |
| **검증** | primary key 중복 0, FK orphan 0, 모든 record의 `accession_class` 채워짐, schema 검증(pandera 권장) 통과 |

```python
# 정규화 핵심: accession_class 도출
import pandas as pd
def accession_class(acc: str) -> str:
    return {"NM": "validated_mrna", "XM": "predicted_mrna",
            "NR": "validated_ncrna", "XR": "predicted_ncrna"}.get(acc[:2], "other")
```

---

### Stage 4: Orthology backbone from `gene_orthologs.gz`

| 항목 | 내용 |
|---|---|
| **목표** | 인간 anchor로 14 target species에 대한 pairwise ortholog 후보 추출 |
| **인풋** | `raw/ncbi_bulk/gene_orthologs.*.gz`, `species_manifest` (14 taxid) |
| **도구** | Python (pandas) |
| **액션** | gene_orthologs.gz 파싱 → `(tax_id=9606, GeneID=h, Other_tax_id ∈ {13 target taxids}, Other_GeneID=t)` 행만 추출. 반대 방향(`tax_id ∈ targets, Other_tax_id=9606`)도 수집해 양방향 edge 구성 |
| **아웃풋** | `orthology/ortholog_candidates.parquet`: 컬럼 `(human_GeneID, target_taxid, target_GeneID, direction)` |
| **검증** | 인간 anchor 컬럼이 모두 9606, target_taxid가 14 species 중 하나로 한정, edge 중복 제거 후 row count 로깅 |

```python
import gzip, pandas as pd
TARGET_TAXIDS = {9598,9597,9595,9601,9600,9590,9541,9545,9483,27679,9470,9447,110931}
rows = []
with gzip.open("raw/ncbi_bulk/gene_orthologs.<date>.gz", "rt") as fh:
    next(fh)  # header
    for line in fh:
        tax_id, GeneID, rel, Other_tax_id, Other_GeneID = line.rstrip("\n").split("\t")
        if rel != "Ortholog": continue
        if int(tax_id) == 9606 and int(Other_tax_id) in TARGET_TAXIDS:
            rows.append((int(GeneID), int(Other_tax_id), int(Other_GeneID), "h_to_t"))
        elif int(Other_tax_id) == 9606 and int(tax_id) in TARGET_TAXIDS:
            rows.append((int(Other_GeneID), int(tax_id), int(GeneID), "t_to_h"))
df = pd.DataFrame(rows, columns=["human_GeneID","target_taxid","target_GeneID","direction"])
df.drop_duplicates(["human_GeneID","target_taxid","target_GeneID"]).to_parquet(
    "orthology/ortholog_candidates.parquet")
```

---

### Stage 5: Strict singleton + human paralog 배제 + MT 배제

| 항목 | 내용 |
|---|---|
| **목표** | (i) MT/non-protein-coding 배제, (ii) human paralog 발견되는 모든 family 배제(many:1), (iii) 14 species coverage + per-species 1 GeneID 검증 |
| **인풋** | `ortholog_candidates.parquet`, `gene_index.parquet` (chromosome, gene_type) |
| **도구** | Python (pandas + 그래프 분석) |
| **액션** | **Step A (MT/biotype 사전 필터)**: 인간 GeneID 중 `chromosome == 'MT'` 또는 `type_of_gene != 'protein-coding'`인 anchor 제거. 또한 어느 target species에서든 매핑된 target GeneID가 MT/non-coding이면 그 edge 제거.<br>**Step B (many:1 = human paralog 배제)**: 각 `(target_taxid, target_GeneID)`에 연결된 distinct `human_GeneID` 수가 2 이상이면 → 거기 연결된 모든 human anchor를 `human_paralog_detected`로 마킹해 배제.<br>**Step C (1:many forward 배제)**: 각 `(human_GeneID, target_taxid)`에 distinct `target_GeneID`가 2 이상이면 → `multi_gene_target_species`로 배제.<br>**Step D (14-species coverage)**: 남은 human anchor에 대해 target_taxid의 distinct count == 13(인간 제외)이어야 통과. 부족하면 `missing_taxon` |
| **아웃풋** | `orthology/strict_singleton.parquet` (통과 set: human_GeneID + 13 target GeneIDs), `orthology/rejected.parquet` (사유 컬럼: `mt_gene`/`non_coding`/`human_paralog_detected`/`multi_gene_target_species`/`missing_taxon`) |
| **검증** | 통과 set 모두 14 species × 1 GeneID, rejection 사유 100% 분류, lineage-expansion family는 Step B/C에서 자동 reject됨을 확인 |

```python
# Step B 핵심 (many:1 paralog 검출)
target_human_counts = (df.groupby(["target_taxid","target_GeneID"])
                         ["human_GeneID"].nunique())
many_to_one_targets = target_human_counts[target_human_counts > 1].index
paralog_humans = set(df.merge(
    pd.DataFrame(list(many_to_one_targets), columns=["target_taxid","target_GeneID"]),
    on=["target_taxid","target_GeneID"])["human_GeneID"])
df_clean = df[~df["human_GeneID"].isin(paralog_humans)].copy()
```

**핵심 통찰**: olfactory receptor, MHC, KIR 같은 lineage-expanded family는 거의 모두 Step B 또는 Step C에서 자동 reject되므로 별도 블랙리스트 불필요(정책 #6).

---

### Stage 6: 종별 대표 CDS 선택 (NM_ 우선, XM_ fallback)

| 항목 | 내용 |
|---|---|
| **목표** | 통과한 모든 (human anchor, species) 조합에 대해 codon analysis용 대표 transcript/CDS 1개씩 선택 |
| **인풋** | `strict_singleton.parquet`, `product_index.parquet` |
| **도구** | Python |
| **액션** | 종별로 다음 결정 트리:<br>① `accession_class == validated_mrna` (NM_) 후보 존재 시 NM_만 사용<br>② 없으면 XM_ 후보로 진행<br>③ `transcript_select_category` (MANE/RefSeq Select) 가용 시 우선<br>④ `protein_length` 최대 우선<br>⑤ `transcript_length` 최대 우선<br>⑥ `transcript_accession` lexicographic 최소 (deterministic tie-break)<br>단, `partial=true` 또는 `cds_accession` 결손 후보는 제외 |
| **아웃풋** | `selection/representative_cds.parquet` (human_GeneID, taxid, GeneID, transcript_accession, cds_accession, protein_accession, protein_length, accession_class, selection_rule_id), `selection_audit.parquet` |
| **검증** | family×species마다 정확히 1 record. NM_ 우선 적용 확인 (species별 NM_/XM_ 비율 통계 출력) |

---

### Stage 7: CDS 무결성 QC

| 항목 | 내용 |
|---|---|
| **목표** | alignment 입력 적합성 확인. PAGAN2가 frame correction 안 하므로 사전 검사 엄격 |
| **인풋** | `representative_cds.parquet` + `cds.fna` 파일들 |
| **도구** | Python (Biopython) |
| **액션** | 각 CDS에 대해:<br>(1) 길이 ≥ 60 nt<br>(2) 길이 % 3 == 0<br>(3) 시작 코돈 ATG 또는 alt(GTG/TTG/CTG; 옵션)<br>(4) terminal stop 존재 시 alignment 입력에서 **제거하고** metadata에 `had_terminal_stop=true` 기록 (정책 #7 trimming 없음과 별개; terminal stop 제거는 aligner 호환 목적)<br>(5) **internal stop = 0** (검출 시 family fail)<br>(6) N 비율 ≤ 1%<br>(7) IUPAC 외 문자 0 |
| **아웃풋** | `qc/cds_qc.parquet` (passed/failed + 사유), 통과 CDS만 다음 stage로 |
| **검증** | family 단위로 14/14 CDS가 통과해야 그 family가 살아남음. 한 species라도 fail → family 통째 reject(`no_valid_cds`) |

```python
from Bio.Data.CodonTable import unambiguous_dna_by_id
STOP_CODONS = set(unambiguous_dna_by_id[1].stop_codons)  # standard table
def has_internal_stop(seq: str) -> bool:
    for i in range(0, len(seq)-3, 3):  # exclude terminal codon
        if seq[i:i+3] in STOP_CODONS:
            return True
    return False
```

---

### Stage 8: 종간 sanity check

| 항목 | 내용 |
|---|---|
| **목표** | family 내 14 CDS가 서로 합리적인 범위에 있는지 확인 (mis-orthology/paralog 잔존 마지막 방어선) |
| **인풋** | Stage 7 통과 CDS |
| **도구** | Python |
| **액션** | family별로 (a) protein_length의 max/min ratio, (b) 빠른 pairwise protein identity sketch (k-mer 또는 짧은 mafft `--auto` 후 columnwise identity 평균) |
| **아웃풋** | `qc/family_sanity.parquet` — ratio>1.5는 `flag`, ratio>2.0은 `fail`. identity median<40%는 `flag` |
| **검증** | flag된 family는 통과하되 alignment 후 별도 검토 권장. fail family는 reject |

---

### Stage 9: SYMBOL.fasta 생성

| 항목 | 내용 |
|---|---|
| **목표** | 통과한 family별로 `fastas/{HUMAN_SYMBOL}.fasta` (14 sequence, header = species token) 생성 |
| **인풋** | Stage 7/8 통과 family + representative CDS sequence |
| **도구** | Python (Biopython) |
| **액션** | (1) human symbol 추출(gene_index의 `symbol`), 파일 시스템 비안전 문자(`/`, `:`, 공백 등)를 `_`로 치환하되 원본 symbol은 meta.tsv에 보존<br>(2) terminal stop 제거된 CDS sequence 사용<br>(3) header는 species token만 (`>human`, `>chimpanzee`, ...) — 공백 없음<br>(4) meta.tsv에 taxid/GeneID/transcript_accession/cds_accession/protein_length/accession_class/had_terminal_stop 기록 |
| **아웃풋** | `fastas/{SYMBOL}.fasta`, `fastas/{SYMBOL}.meta.tsv` |
| **검증** | 각 FASTA의 sequence 개수 == 14, header set == 14 token, filename 충돌 0 |

**FASTA 예시**:
```
>human
ATGGAT...
>chimpanzee
ATGGAT...
>bonobo
ATGGAT...
...
>Philippine_flying_lemur
ATGGAT...
```

**meta.tsv 필수 컬럼**:
- `human_symbol`, `human_symbol_original`, `human_GeneID`
- `taxid`, `token`, `scientific_name`
- `GeneID`, `transcript_accession`, `cds_accession`, `protein_accession`
- `protein_length`, `cds_length`, `accession_class`, `selection_rule_id`
- `had_terminal_stop`

---

### Stage 10: ROADIES tree 검증 (사용자 제공)

| 항목 | 내용 |
|---|---|
| **목표** | 사용자가 제공한 rooted binary tree가 PAGAN2 입력 요구사항을 충족하는지 검증만 수행 (rooting/binarization 로직은 파이프라인에 없음 — 정책 #3) |
| **인풋** | `trees/species_tree.rooted.binary.nwk` |
| **도구** | dendropy |
| **액션** | (1) leaf set == 14 species token (Philippine_flying_lemur 포함)<br>(2) `is_rooted == true`<br>(3) 모든 internal node가 정확히 2개 child<br>(4) branch length 0이거나 결손인 가지 → PAGAN2 호환을 위해 작은 epsilon(예: 1e-6) 대체본 별도 저장<br>(5) flying lemur가 outgroup 위치인지 확인 (가장 깊은 split이 flying lemur vs rest 13 primates인가) |
| **아웃풋** | 검증 통과 시 tree 그대로 사용. epsilon 대체본은 `species_tree.rooted.binary.pagan_safe.nwk`로 저장 |
| **검증** | 5개 검사 모두 통과해야 Stage 11 진입 |

---

### Stage 11: PAGAN2 codon-aware alignment

| 항목 | 내용 |
|---|---|
| **목표** | family별 codon-aware alignment 생성. PAGAN2 primary + PRANK fallback. 정책 #7에 따라 untrimmed만 |
| **인풋** | `fastas/{SYMBOL}.fasta`, `trees/species_tree.rooted.binary.pagan_safe.nwk` |
| **도구** | PAGAN2 (Apptainer), PRANK (conda), MACSE (옵션) |
| **액션** | Stage 7이 frame issue를 미리 잘라냈으므로 거의 모든 family가 PAGAN2로 직진. 라우팅:<br>(1) PAGAN2 실행 → 성공 + Stage 12 QC 통과 → 채택<br>(2) PAGAN2 실패 또는 QC 실패 → PRANK 재실행 → 성공 + QC 통과 → 채택<br>(3) 둘 다 실패 → `alignment_failed` 사유와 함께 reject |
| **아웃풋** | `alignments/pagan2/{SYMBOL}.fas` 또는 `alignments/prank/{SYMBOL}.best.fas`, `alignments/aligner_log.parquet` (어느 aligner가 어느 family에 성공했는지) |
| **검증** | alignment 파일 존재, sequence 개수 14, length % 3 == 0 |

```bash
# PAGAN2 (Apptainer wrapper)
bin/pagan2-run \
  --seqfile /data/fastas/${SYMBOL}.fasta \
  --treefile /data/trees/species_tree.rooted.binary.pagan_safe.nwk \
  --codons \
  --outfile /data/alignments/pagan2/${SYMBOL} \
  --outformat fasta \
  --config-log-file /data/alignments/pagan2/${SYMBOL}.cfg

# PRANK fallback
prank \
  -d=fastas/${SYMBOL}.fasta \
  -t=trees/species_tree.rooted.binary.pagan_safe.nwk \
  -codon -prunedata -prunetree \
  -o=alignments/prank/${SYMBOL}
```

**PAGAN2 주의사항**:
- codon mode는 **first reading frame 가정 + frame correction 없음** → Stage 7 QC가 필수
- terminal stop은 Stage 7에서 이미 제거됨 → aligner 비교 일관성 확보
- guide tree는 rooted binary 필수 → Stage 10이 보장

---

### Stage 12: Alignment QC

| 항목 | 내용 |
|---|---|
| **목표** | 각 alignment의 결과 무결성 확인 |
| **인풋** | `alignments/{pagan2,prank}/*.fas` |
| **도구** | Python (Biopython) |
| **액션** | (1) sequence 개수 == 14<br>(2) alignment length % 3 == 0<br>(3) 코돈 단위로 internal stop 검사 — 발견 시 alignment 폐기 후 PRANK fallback (이미 PRANK였다면 reject)<br>(4) per-sequence gap fraction 산출 (> 50% flag, > 80% fail)<br>(5) header set ≡ tree leaf set 재검증 |
| **아웃풋** | `qc/alignment_qc.parquet` (status: `pass`/`flag`/`fail`) |
| **검증** | family당 status 1개, fail 사유 100% 카테고리화 |

---

### Stage 13: 최종 리포트 + provenance

| 항목 | 내용 |
|---|---|
| **목표** | 입력 버전·NCBI snapshot·rejection 통계·최종 family 수를 단일 HTML로 요약 |
| **인풋** | 모든 stage의 산출물 |
| **도구** | Python (pandas + jinja2) |
| **액션** | (1) NCBI snapshot date·14 GCF 버전·CLI 버전 명시<br>(2) stage별 입출력 row count<br>(3) rejection 사유 분포 그래프<br>(4) 종별 NM_/XM_ 비율<br>(5) X-linked vs autosomal singleton 비율 (Y는 별도 표기)<br>(6) PAGAN2/PRANK 성공률 |
| **아웃풋** | `reports/summary.html`, `reports/provenance.html` |
| **검증** | 모든 통계가 raw 산출물에서 자동 재계산되는지 |

---

## 4. Snakemake workflow 스켈레톤

```python
# workflow/Snakefile
configfile: "config/pipeline_config.yaml"
SPECIES = config["species"]           # list of tokens
TAXIDS = config["taxids"]              # list of 14 taxids

rule all:
    input: "reports/summary.html"

rule preflight:
    input:
        manifest="config/species_manifest.tsv",
        tree="trees/species_tree.rooted.binary.nwk"
    output: "reports/preflight.json"
    conda: "envs/ncbi.yml"
    script: "scripts/00_preflight.py"

rule bulk_freeze:
    output:
        orth="raw/ncbi_bulk/gene_orthologs.gz",
        info="raw/ncbi_bulk/gene_info.gz",
        chk="raw/ncbi_bulk/MANIFEST.checksums"
    shell: "scripts/01_bulk_freeze.sh"

rule gene_package:
    input: rules.preflight.output
    output: directory("raw/gene_packages/{taxid}")
    conda: "envs/ncbi.yml"
    shell: "scripts/02_gene_package.sh {wildcards.taxid} {output}"

rule normalize_indexes:
    input:
        packages=expand("raw/gene_packages/{taxid}", taxid=TAXIDS),
        info=rules.bulk_freeze.output.info
    output:
        gene="indexes/gene_index.parquet",
        product="indexes/product_index.parquet",
        cds="indexes/cds_index.parquet"
    script: "scripts/03_normalize.py"

rule orthology_backbone:
    input:
        orth=rules.bulk_freeze.output.orth,
        gene=rules.normalize_indexes.output.gene
    output: "orthology/ortholog_candidates.parquet"
    script: "scripts/04_backbone.py"

rule strict_singleton:
    input:
        cands=rules.orthology_backbone.output,
        gene=rules.normalize_indexes.output.gene
    output:
        passed="orthology/strict_singleton.parquet",
        rejected="orthology/rejected.parquet"
    script: "scripts/05_singleton.py"  # MT 배제 + many:1 paralog 배제 + 14-coverage

rule representative_cds:
    input:
        singleton=rules.strict_singleton.output.passed,
        product=rules.normalize_indexes.output.product
    output: "selection/representative_cds.parquet"
    script: "scripts/06_select_cds.py"  # NM_ 우선, XM_ fallback

rule cds_qc:
    input:
        sel=rules.representative_cds.output,
        cds=rules.normalize_indexes.output.cds
    output: "qc/cds_qc.parquet"
    script: "scripts/07_cds_qc.py"

rule family_sanity:
    input: rules.cds_qc.output
    output: "qc/family_sanity.parquet"
    script: "scripts/08_family_sanity.py"

rule write_fastas:
    input:
        qc=rules.cds_qc.output,
        sanity=rules.family_sanity.output,
        sel=rules.representative_cds.output
    output: directory("fastas")
    script: "scripts/09_write_fastas.py"

rule tree_validate:
    input: "trees/species_tree.rooted.binary.nwk"
    output: "trees/species_tree.rooted.binary.pagan_safe.nwk"
    script: "scripts/10_tree_validate.py"

def family_targets(_):
    import glob
    return [f"alignments/final/{p}.fas"
            for p in (s.split("/")[-1].replace(".fasta","")
                      for s in glob.glob("fastas/*.fasta"))]

rule align_pagan2:
    input:
        fa="fastas/{symbol}.fasta",
        tree=rules.tree_validate.output
    output: "alignments/pagan2/{symbol}.fas"
    log: "logs/pagan2/{symbol}.log"
    shell:
        r"""
        mkdir -p alignments/pagan2 logs/pagan2
        bin/pagan2-run \
          --seqfile /data/{input.fa} \
          --treefile /data/{input.tree} \
          --codons --outfile /data/alignments/pagan2/{wildcards.symbol} \
          --outformat fasta > {log} 2>&1 || true
        test -s {output}
        """

rule align_prank_fallback:
    input:
        fa="fastas/{symbol}.fasta",
        tree=rules.tree_validate.output
    output: "alignments/prank/{symbol}.best.fas"
    log: "logs/prank/{symbol}.log"
    conda: "envs/alignment.yml"
    shell:
        r"""
        mkdir -p alignments/prank logs/prank
        prank -d={input.fa} -t={input.tree} -codon -prunedata -prunetree \
              -o=alignments/prank/{wildcards.symbol} > {log} 2>&1
        """

rule align_final:
    input: pagan="alignments/pagan2/{symbol}.fas"
    output: "alignments/final/{symbol}.fas"
    script: "scripts/11_route_align.py"  # PAGAN2 결과 QC → 실패 시 PRANK 트리거

rule alignment_qc:
    input: family_targets
    output: "qc/alignment_qc.parquet"
    script: "scripts/12_aln_qc.py"

rule summary:
    input: rules.alignment_qc.output
    output: "reports/summary.html"
    script: "scripts/13_report.py"
```

---

## 5. 리스크 매트릭스

| 리스크 | 대응 |
|---|---|
| GCF 버전 ≠ gene package "current" annotation release | Stage 0에서 release_id 차이 자동 감지, 차이 큰 species는 genome+GFF 경로로 fallback |
| 인간 T2T-CHM13v2.0의 transcript accession이 GRCh38과 다름 | Stage 6에서 species별 NM_/XM_ 비율 로깅. 인간 NM_가 비정상적으로 적으면 release 확인 |
| 일부 GCF가 Y chromosome 미포함 | Stage 0의 `has_Y` 표로 사전 가시화. Y-linked singleton 0은 자연스러운 결과로 reporting |
| `gene_orthologs.gz` snapshot 시점과 gene package 시점 불일치 | 두 다운로드를 같은 wall-clock 시간대(수 시간 이내)에 수행, snapshot date 함께 기록 |
| PAGAN2 Docker pull 실패/HPC에서 Docker 금지 | Apptainer wrapper로 우회. `bin/pagan2-run`이 두 backend 모두 추상화 |
| PAGAN2가 frame correction 안 함 | Stage 7 QC를 통과한 CDS만 진입 → 알고리즘적 보호 |
| Macaca fascicularis GCF의 `-RS_YYYY_MM` annotation suffix 형식 | manifest에는 `GCF_037993035.2`로만 기록. Stage 0에서 release_id 별도 컬럼 |
| 사용자 제공 tree가 unrooted/multifurcating으로 잘못 제공됨 | Stage 10에서 즉시 fail. 사용자에게 명확한 에러 메시지 |

---

## 6. 즉시 시작할 수 있는 명령어 (3시간 내 MVP 0 도달)

```bash
# 1) repo skeleton
mkdir -p config raw/ncbi_bulk raw/gene_packages indexes orthology selection \
         qc fastas trees alignments/pagan2 reports workflow/scripts \
         containers bin envs logs

# 2) conda envs
conda env create -f envs/ncbi.yml
conda env create -f envs/alignment.yml

# 3) PAGAN2 컨테이너 (Docker 또는 Apptainer 중 가능한 쪽)
apptainer build containers/pagan2.sif docker://ariloytynoja/pagan2

# 4) species_manifest.tsv 작성 (위 14종 표 그대로)

# 5) Pre-flight 수동 실행 (Snakemake 작성 전 sanity check)
conda activate ncbi
bash workflow/scripts/01_bulk_freeze.sh         # gene_orthologs.gz + gene_info.gz 다운로드
python workflow/scripts/00_preflight.py         # 14 GCF release_id, has_Y, tree validation
```

이 시점에 `reports/preflight.json`에서 14/14 검사 통과를 보면 Stage 1~13을 Snakemake로 일괄 실행할 수 있는 상태가 됩니다.

---

## 7. 다음 작업 권장 순서

| 마일스톤 | 범위 | 완료 기준 |
|---|---|---|
| **MVP 0** | §6의 6개 명령어 실행 | `reports/preflight.json`에서 14/14 검사 통과 |
| **MVP 1** | Stage 1~3 스크립트 작성 | 14종 gene package + 정규화 index 생성 |
| **MVP 2** | Stage 4~5 스크립트 작성 | `strict_singleton.parquet` 생성, 사유 분포 확인. 이 시점에 "최종 family 수"의 대략적 추정 가능 |
| **MVP 3** | Stage 6~9 스크립트 작성 | SYMBOL.fasta 생성. 14종 1:1 set 첫 산출물 손에 들어옴 |
| **V1** | Stage 10~13 스크립트 작성 | ROADIES tree 검증 + PAGAN2 정렬 + 최종 리포트 |

각 MVP마다 reject 사유 분포와 family 수를 확인하면, downstream 분석(selection, gene tree, ASR)에 충분한 수의 family가 확보되는지 조기에 판단할 수 있습니다.

---

## 부록 A: 핵심 KPI

| 카테고리 | 지표 | 목표 |
|---|---|---|
| 정합성 | strict singleton false positive | 0 |
| 정합성 | 종별 대표 CDS 수 | family×species 모두 정확히 1 |
| 정합성 | FASTA header / tree leaf mismatch | 0 |
| 정합성 | provenance 추적 완전성 | 100% (taxid, GeneID, accession, GCF, snapshot date 모두) |
| 정합성 | 동일 input·snapshot에서 재현성 | bit-identical manifest |
| 운영 | NCBI download 성공률 (with retry) | 100% within window |
| 운영 | Schema validation 실패 | 0 |
| 운영 | Citation manifest 완전성 | 100% |

## 부록 B: rejection 사유 카테고리 (완전 분류)

| 사유 | 발생 stage | 의미 |
|---|---|---|
| `mt_gene` | Stage 5 | 인간 또는 어느 target species에서든 MT 유전자 |
| `non_coding` | Stage 5 | `type_of_gene != 'protein-coding'` |
| `human_paralog_detected` | Stage 5 | many:1 (target GeneID가 여러 human GeneID와 ortholog) |
| `multi_gene_target_species` | Stage 5 | 1:many (한 human GeneID가 target species에서 여러 GeneID 매핑) |
| `missing_taxon` | Stage 5 | 14 target species 중 1개 이상 누락 |
| `no_valid_transcript` | Stage 6 | NM_, XM_ 모두 부재 또는 모두 partial=true |
| `cds_qc_fail` | Stage 7 | 길이/start/stop/N/internal stop 위반 |
| `family_length_outlier` | Stage 8 | length max/min ratio > 2.0 |
| `family_identity_low` | Stage 8 | median pairwise protein identity < 임계 |
| `alignment_failed` | Stage 11 | PAGAN2와 PRANK 모두 실패 |
| `alignment_qc_fail` | Stage 12 | length/stop/gap 위반 |

---

**참고 문서**:
- NCBI Datasets CLI: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/
- NCBI Gene FTP: https://ftp.ncbi.nlm.nih.gov/gene/DATA/
- PAGAN2: https://github.com/ariloytynoja/pagan-msa
- PRANK: http://wasabiapp.org/software/prank/
- MACSE: https://www.agap-ge2pop.org/macse/
- ROADIES: https://turakhia.ucsd.edu/ROADIES/
