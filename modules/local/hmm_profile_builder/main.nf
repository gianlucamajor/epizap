
process HMM_PROFILE_BUILDER {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    label "few_cpu"
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'ignore' } 
    maxRetries 3

    input:
    path(msaFile)
    
    output:
    path("hmm/*.hmm")

    script:
    """
    mkdir hmm
    hmmbuild hmm/${msaFile.baseName}.hmm ${msaFile}
    """
}
