process GRAPH_UPDATER {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${metaGraph.id}"
    label "few_cpu"
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    tuple val(metaGraph), path(graphFile)
    tuple val(metaReport), path(reportFile)
    tuple val(metaFeature), path(featureFile)

    output:
    path("graph-updated/*.pickle"), emit: graph


    script:
    """
    mkdir graph-updated 
    graph_updater.py ${graphFile} ${reportFile} -f ${featureFile} -o graph-updated
    """
    

}