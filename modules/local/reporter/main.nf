process EPITOPE_REPORTER {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${meta.id}"
    label "med_cpu"
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    tuple val(meta), path(graphFile)

    output:
    path("epitopes_report/*.json"), emit: jsonReport
    path("epitopes_report/*.fasta"), emit: fastaReport
    path("epitopes_report/*.pickle"), emit: graphReport

    script:
    """
    mkdir -p epitopes_report
    epitope_reporter.py ${graphFile} -o epitopes_report/
    cp ${graphFile} epitopes_report/

    """
}
// python3 bin/epitope_reporter.py results/graph-updated/CCC_severe_a_mapped-segment-graph-graph-cc-id-msa-epitopes.pickle -o results/ept-reported/