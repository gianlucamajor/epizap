process LONELY {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${file.name}"
    label "one_cpu"
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'ignore' } 
    maxRetries 3

    input:
    tuple val(meta), path(file)

    output:
    path("lonely/*.lone")

    script:
    """
    mkdir lonely
    cp  ${file} lonely/${file.baseName}.lone
    """
}