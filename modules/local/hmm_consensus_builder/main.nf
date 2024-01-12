process HMM_CONSENSUS_BUILDER {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    label "one_cpu"
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    path(hmmFile)

    output:
    path("consensus/*.consensus")


    script:
    """
    mkdir consensus
    hmmemit -c ${hmmFile} > consensus/${hmmFile.baseName}.consensus
    """
}