
process HMM_PROFILE_BUILDER {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${msaFile.baseName}"
    label "few_cpu"
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    path(msaFile)
    
    output:
    path("hmm/*.hmm")

    script:
    """
    mkdir hmm
    hmmbuild --cpu $task.cpus hmm/${msaFile.baseName}.hmm ${msaFile}
    """
}
