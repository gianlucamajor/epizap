process IDX_BUILDER {
    publishDir "${params.outdir}", mode: 'copy', overwrite: false
    tag "${meta.id}"
    cpus 6
    input:
    tuple val(meta), path(fastaRef)

    output:
    tuple val(meta), path("${fastaRef.simpleName}-idx"), emit: idx

    script:
    """
    mkdir "${fastaRef.simpleName}-idx"
    bowtie2-build --threads $task.cpus ${fastaRef} "${fastaRef.simpleName}-idx/${fastaRef.simpleName}"
    """
}