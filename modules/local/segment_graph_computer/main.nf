process SEGMENT_GRAPH_COMPUTER {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${metaSeg.id} - threshold ${params.graph_segment_threshold}"
    label "med_cpu_high_memory"
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    tuple val(metaSeg), path(segmentFile)

    output:
    path("graph/*.pickle"), emit: segmentGraph

    script:
    """
    mkdir graph
    segment_graph_computer.py ${segmentFile} -t ${params.graph_segment_threshold} -o graph -p ${task.cpus}
    """
}