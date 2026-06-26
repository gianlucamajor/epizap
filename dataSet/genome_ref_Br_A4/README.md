### Reference Genome Dataset

This directory contains the reference genome sequence and annotation files for the *Trypanosoma cruzi* strain **Brazil Clone A4**. These reference files are essential for genomic mapping, alignment, and downstream bioinformatics analyses.

**Assembly Details:**
* **Organism:** *Trypanosoma cruzi*
* **Strain:** Brazil Clone A4
* **Assembly Name:** ASM1503362v1
* **NCBI Accession:** [GCA_015033625.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_015033625.1/)

---

### Files & Data Integrity

The table below provides the details of the available compressed files, including their descriptions and MD5 checksums for data integrity verification before extraction and usage.

| File Name | MD5 Checksum | Description |
| :--- | :--- | :--- |
| `t_cruzi_br_a4.fna.zip` | `427cfd0bc3c3728723c8a6b9a741de4d` | Genomic FASTA format (`.fna`) file containing the complete nucleotide sequence of the assembly (compressed). |
| `genomic.gff.zip` | `22dd7ce4d8cd52a3d1a1af5b58148f41` | General Feature Format (`.gff`) file containing the genomic annotations and feature coordinates (compressed). |

---

### Usage Note
Please ensure you unzip the `.zip` files to access the raw `.fna` and `.gff` formats required by most bioinformatics tools (such as BLAST, BWA, or bedtools) prior to running your pipeline.