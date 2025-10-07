# epizap

Epitope predictor

### Requirements

#### NEXTFLOW

-> Nextflow can be used on any POSIX-compatible system (Linux, macOS, etc), and on Windows through WSL. It requires Bash 3.2 (or later) and Java 11 (or later, up to 23). You can see which version you have using the following command:
java -version

More details on: https://www.nextflow.io/docs/latest/install.html

#### DOCKER

https://docs.docker.com/get-started/

#### Image from docker Hub

https://hub.docker.com/repository/docker/gianlucamajor/epizap/general

##### Current version

docker pull gianlucamajor/epizap:v0.2.0-alpha

### Quick Start 


```
nextflow epizap.nf --reads dataSet/k-dna-insert-reads/CCC_severe_a.fastq
```

#### Expected result
```
N E X T F L O W   ~  version 24.04.4

Launching `epizap.nf` [tender_mandelbrot] DSL2 - revision: 25b3109834

===================================
Reads: dataSet/k-dna-insert-reads/CCC_severe_a.fastq
Peptides: /home/gianluca/workspace/epizap/dataSet/peptides/kDNA-aa-in-frame.fasta
Reference: /home/gianluca/workspace/epizap/dataSet/ref/y_strain_minicircles.fasta
Outdir: results
Peptiede Segment: /home/gianluca/workspace/epizap/results/peptides-segment/*.fasta
Mapper file:

===================================

executor >  local (68)
[de/70bcce] mapper:IDX_BUILDER (y_strain_minicircles)                                  [100%] 1 of 1 ✔
[6b/a92972] mapper:MAPPER (CCC_severe_a)                                               [100%] 1 of 1 ✔
[b7/da805e] mapper:MAX_SCORE_MAPPING_EXTRACTOR (CCC_severe_a)                          [100%] 1 of 1 ✔
[a0/5db3a7] segmentRetriever:SEGMENT_EXTRACTOR (CCC_severe_a)                          [100%] 1 of 1 ✔
[ba/6b3575] segmentRetriever:SEQUENCE_EXTRACTOR (CCC_severe_a)                         [100%] 1 of 1 ✔
[f3/5e3804] msa:MSA (Tc_minicircle_Y_4.93-49-152-34-5.fasta; align)                    [100%] 30 of 30 ✔
[24/e75ed6] mview:MVIEW (Tc-minicircle-Y-4.93-49-152-34-5)                             [100%] 30 of 30 ✔
[6b/7358ae] predicator:PREDICATOR                                                      [100%] 1 of 1 ✔
[e7/b09255] lonelyPeptideRetriever:LONELY (number of segments with lonely peptide: 59) [100%] 1 of 1 ✔
[4b/baa848] msaCoreAndLonelyEpitopesJoiner:JOIN_MSA_CORE_AND_LONELY_EPITOPES           [100%] 1 of 1 ✔

```
---

Other infos 

### mview

https://desmid.github.io/mview/index.html

mview -in fasta -html head -css on -coloring consensus -threshold 100 -consensus on -con_threshold 100 ../msa/control_and_chagasic_patients_WNWZ01000122.1-6914-7101-40_26-msa.msc > WNWZ01000122.1-6914-7101-40_26-msa.html

### IGV-webapp

npx http-server -a localhost /home/gianluca/workspace/igv-webapp/
http://localhost:8080/igv-webapp/?locus=CM026600.1:31354-31391

### Extract features from GFF

bedtools intersect -a /home/gianluca/workspace/epizap/results_21_11_2024/segments/control_and_chagasic_patients_mapped-segment.tsv -b genomic.gff -wb | grep "CDS" | less -S

### Run with mapper file

#nextflow epizap.nf --mapper /home/gianluca/workspace/epizap/results-04-08-2024/only-primary-mapping/control_and_chagasic_patients_opm.bam --ref dataSet/ref/t_cruzi_br_a4.fna  --outdir results-04-08-2024/ --peptides dataSet/peptides/gDNA-kDNA-aa-in-frame-60-ident.fasta --reads /work/data/epizap/gDNA-kDNA-insert-reads/control_and_chagasic_patients.fastq

### Run


nextflow epizap.nf \
    --reads /home/gianluca/workspace/data/epizap/gDNA-kDNA-insert-reads/control_and_chagasic_patients.fastq \
    --peptides /home/gianluca/workspace/data/epizap/peptides/gDNA-kDNA-aa-in-frame-60-ident.fasta \
    --ref /home/gianluca/workspace/data/epizap/ref/t_cruzi_br_a4.fna \
    --annotation /home/gianluca/workspace/data/epizap/ref/genomic.gff \
    --outdir "results_18_09_2025"


### RUN CD-HIT VERIFY EPITOPE FINAL LIST RESULT 

cd-hit -i epitopes-cc-graph.fasta -o epitopes-cc-graph-cdhit-100 -c 1.00 -l 4 -n 5 -d 0 -T 10 -M 2000
#### Given tha the min length of epitope predicted by epizap is 5 than is important define -l parameter as 4. 

-c 1.0, means 100% identity, is the clustering threshold
-l	length of throw_away_sequences, default 10
-n 5 is the word size
-d 0 use sequence name in fasta header till the first white space
-M 16000, to use 16GB RAM
-T 8, to use 8 threads

#### 
Export fastq from bam
samtools fastq control_and_chagasic_patients_msm_sorted.bam -F 256 > all-inserts-mapped-BrA4.fastq
Keep just id on peptides sequence header
seqkit seq -i gDNA-kDNA-aa-in-frame-60-ident.fasta > t-cruzi-all-peptides.fasta

### Epitope reporter 
python3 bin/epitope_reporter.py results_25_09_2025/graph-updated/control_and_chagasic_patients_mapped-segment-graph-graph-cc-id-msa-epitopes.pickle \
    -o results_25_09_2025/ept-reported/ \
    --iedb \
    --iedb-tcruzi-epitopes-hits results_25_09_2025/blastp-tcruzi-iedb-epitopes/EPZ-IEDB.tsv \
    --iedb-tcruzi-epitopes dataSet/iedb/tcruzi/epitope_table_export_1747317719_15052025.csv \
    --iedb-human-epitopes-hits results_25_09_2025/blastp-human-iedb-epitopes/EPZ-HUMAN-IEDB-25-09-25.tsv \
    --iedb-human-epitopes dataSet/iedb/human/human_epitope_table_export_1757095864.csv
    