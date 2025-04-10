process PREDICATOR {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    // tag "${meta.id}"
    label "med_cpu"
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    path(mviewFiles)

    output:
    path("epitopes-by-msa-core/*.fasta"), emit: pepSeg
    path("epitopes-by-msa-core/*.tsv"), emit: msaEpitopeReport

    script:
    """
    mkdir epitopes-by-msa-core
    epitope_predictor.py  ${mviewFiles} -o epitopes-by-msa-core
    """
}