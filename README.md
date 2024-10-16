# epizap
Epitope predictor

### mview
https://desmid.github.io/mview/index.html

mview -in fasta -html head -css on -coloring consensus -threshold 100 -consensus on -con_threshold 100 ../msa/control_and_chagasic_patients_WNWZ01000122.1-6914-7101-40_26-msa.msc > WNWZ01000122.1-6914-7101-40_26-msa.html

### IGV-webapp
npx http-server -a localhost /home/gianluca/workspace/igv-webapp/



### Run with mapper file 
 #nextflow epizap.nf --mapper /home/gianluca/workspace/epizap/results-04-08-2024/only-primary-mapping/control_and_chagasic_patients_opm.bam --ref dataSet/ref/t_cruzi_br_a4.fna  --outdir results-04-08-2024/ --peptides dataSet/peptides/gDNA-kDNA-aa-in-frame-60-ident.fasta --reads /work/data/epizap/gDNA-kDNA-insert-reads/control_and_chagasic_patients.fastq

 ### Run 
 nextflow epizap.nf --reads /work/data/epizap/gDNA-kDNA-insert-reads/control_and_chagasic_patients.fastq --peptides dataSet/peptides/gDNA-kDNA-aa-in-frame-60-ident.fasta --ref dataSet/ref/t_cruzi_br_a4.fna --outdir "results_24_08_2024"

### Requirements
#### NEXTFLOW 
-> Nextflow can be used on any POSIX-compatible system (Linux, macOS, etc), and on Windows through WSL. It requires Bash 3.2 (or later) and Java 11 (or later, up to 23). You can see which version you have using the following command:
 java -version

 More details on: https://www.nextflow.io/docs/latest/install.html

