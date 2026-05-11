# refseq2cds Mini Example

This offline fixture is a tiny RefSeq-like dataset for smoke-testing the
`refseq2cds` graph, CDS QC, FASTA writing, and human coordinate matrix logic
without downloading data from NCBI.

It contains three locked species: `human`, `chimpanzee`, and `gorilla`.

Expected behavior:

| Family | Purpose | Expected |
|---|---|---|
| `PASS1` | 3 species, 1 GeneID each, valid CDS | accepted FASTA and matrix |
| `MISSING1` | missing gorilla ortholog edge | rejected |
| `PARALOG1` | gorilla has two GeneIDs in the component | rejected |
| `STOP1` | strict singleton but human CDS has an internal stop | rejected at CDS QC |
| `NEG1` | valid negative-strand human CDS | accepted FASTA and matrix |

Run:

```bash
refseq2cds test --example mini
```

Or run the offline pipeline directly:

```bash
refseq2cds run \
  --manifest examples/mini/manifest.tsv \
  --input-root examples/mini/inputs \
  --output-root /tmp/refseq2cds_mini \
  --offline \
  --steps all \
  --with-matrices \
  --force
```
