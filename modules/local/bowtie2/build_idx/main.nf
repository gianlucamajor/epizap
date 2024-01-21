process IDX_BUILDER {
    publishDir "${params.outdir}", mode: 'copy', overwrite: false
    tag "${meta.id}"
    label 'many_cpu'
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}
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
