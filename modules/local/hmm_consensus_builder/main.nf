process HMM_CONSENSUS_BUILDER {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    label "few_cpu"

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