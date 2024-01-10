process HMM_CONSENSUS_BUILDER {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    label "one_cpu"

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