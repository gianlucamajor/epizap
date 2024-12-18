process JOIN_MSA_CORE_AND_LONELY_EPITOPES {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    label "med_cpu"
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    path(files)

    output:
    path("epitopes-by-msa-core-and-lonely/epitopes-msa-core-and-lonely.fasta"), emit: msaCoreAndLonelyEpitopes


    script:
    """
    mkdir epitopes-by-msa-core-and-lonely
    cat ${files} > epitopes-by-msa-core-and-lonely/epitopes-msa-core-and-lonely.fasta
    """ 
}