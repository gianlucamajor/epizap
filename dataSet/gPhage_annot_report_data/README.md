# gPhage Annotation Report Data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20941568.svg)](https://doi.org/10.5281/zenodo.20941568)

Inputs for the `epizap-annot-report.nf` workflow (`EPITOPE_ANNOT_REPORTER` / `bin/epitope_reporter.py`): BLAST hit tables and IEDB epitope metadata used to annotate phage display-derived epitope candidates.

The `.zip` files are not stored in this repo — download them from Zenodo by running:

```bash
bin/download_gPhage_annot_report_data.sh
```

This downloads each file into this directory, verifies its MD5 checksum, and unzips it.

| File | Description | MD5 |
|---|---|---|
| `proteome_hits.tsv.zip` | BLAST hits of epitopes vs. the *T. cruzi* proteome — maps each candidate to its source protein. | `c2bda4d2d3c17b4cad746259b6d4edf8` |
| `iedb_tcruzi_epitopes_hits.tsv.zip` | BLAST hits vs. the IEDB *T. cruzi* epitope DB — flags matches to known epitopes (paired with the CSV below via `sseqid`). | `e0f710a7d02dfb1602e40b8ecfecad9f` |
| `iedb_tcruzi_epitopes_1747317719_15052025.csv.zip` | IEDB *T. cruzi* epitope metadata export. | `ac09ac6e3e2cacb467ae05e46767ae97` |
| `iedb_human_epitopes_hits.tsv.zip` | BLAST hits vs. the IEDB human epitope DB — flags cross-reactivity with human epitopes (optional). | `1134db0998baa6aea51d927e92f6bf9c` |
| `iedb_human_epitopes_1757095864.csv.zip` | IEDB human epitope metadata export. | `436cffd4c67c5144a8eb4fb61586adb0` |
