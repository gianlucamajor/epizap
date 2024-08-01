process PRE_EPITOPE_RETRIEVER {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    label "med_cpu"
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    path(files)

    output:
    path("pre-epitopes/pre-epitopes.fasta"), emit: preEpitopes


    script:
    """
    mkdir pre-epitopes
    cat ${files} > pre-epitopes/pre-epitopes.fasta
    """ 
}