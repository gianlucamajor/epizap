### Data Availability

All of these pre-processed files and the complete dataset are publicly available for download on Zenodo. You can access the repository via the DOI link below:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20854186.svg)](https://doi.org/10.5281/zenodo.20854186)

### Input Datasets

This directory contains the pre-processed files from the gPhage project, which serve as the input data for the EPIZAP pipeline to predict *Trypanosoma cruzi* epitopes recognized by the antibody response of Chagas disease patients. The dataset comprises eight individual sample files, representing two biological replicates (A and B) across four distinct clinical groups:

* **Control:** Healthy volunteers (Negative for Chagas disease)
* **Asymptomatic:** Patients with no EKG abnormalities
* **Mild CCC:** Patients with mild Chagas cardiomyopathy
* **Severe CCC:** Patients with severe Chagas cardiomyopathy

**List of Individual Sample Files:**
`control_a`, `control_b`, `asympto_a`, `asympto_b`, `CCC_mild_a`, `CCC_mild_b`, `CCC_severe_a`, `CCC_severe_b`

---

### Consolidated Files & File Integrity

To streamline the pipeline execution, data across all patient groups have been merged into consolidated files. The table below provides the file details, total record counts, and MD5 checksums for data integrity verification prior to running EPIZAP.

| File Name | Records | MD5 Checksum | Description |
| :--- | :--- | :--- | :--- |
| `DNA_inserts_all_patients_groups.fastq.zip` | 2,297,587 | `9fcc12426e2ebf6c520edc5dd42e818f` | Contains the merged coding DNA inserts from all eight clinical samples (compressed). |
| `DNA_insert_IDs_by_group.tsv.zip` | 2,297,587 | `785fb95877242e021d3b49772b5209aa` | Specifies which patient clinical group each DNA insert ID belongs to (compressed). |
| `peptides_encoded_by_DNA_inserts.fasta.zip` | 62,282 | `a3ac336dc898e18a5571e5851e40f9bb` | Contains the amino acid sequences of the unique peptides encoded by the DNA inserts (compressed). |
| `peptides_encoded_by_DNA_inserts.tsv.zip` | 62,282 | `29e64a654489eea63a33690d9fe067fe` | A tabular (.tsv) version of the peptide data, including sequence length information (compressed). |

---


### DNA Insert to Peptide Mapping

Importantly, the uniqueness of the peptides within the `peptides_encoded_by_DNA_inserts` is defined by their genomic origin (the specific coding DNA insert) rather than solely by their amino acid sequence. Due to the degeneracy of the genetic code, distinct nucleotide sequences can frequently encode structurally identical amino acid sequences. Therefore, each peptide is assigned a unique identifier (ID) that links it back to its originating DNA insert, preserving full traceability to the corresponding genomic sequence. This strict distinction is maintained throughout the entire processing pipeline.

In practice, this relationship is established by the numeric suffix appended at the end of the DNA insert ID within the FASTQ file.

For example, a DNA insert with the ID `@M01677:61:000000000-AHLMU:1:1106:17399:9980_82653` encodes the peptide identified as `82653` (indicated by the `_82653` suffix). You can find this exact corresponding sequence in the `peptides_encoded_by_DNA_inserts.fasta` file, structured as follows:

```fasta
>82653
TPTDEEKDWRRISRGGQFLAVGRGSTANWSKKPLKQVDQNVGN
