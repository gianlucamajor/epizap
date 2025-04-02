process PEPTIDE_CLUSTERING_CC {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${metaGraph.id}"
    label "few_cpu"
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    tuple val(metaGraph), path(graphFile)
    tuple val(metaPep), path(pepFile)

    output:
    path("peptides-clustering-cc/*.pickle"), emit: segmentGraph
    path("peptides-clustering-cc/connected-components/*.fasta"), emit: peptidesByCC
    path("peptides-clustering-cc/*.fasta"), emit: onlyOnePeptideList

    script:
    """
    mkdir peptides-clustering-cc
    peptide_clustering_by_CC.py ${graphFile} ${pepFile} -o peptides-clustering-cc -cc-output connected-components
    """
}