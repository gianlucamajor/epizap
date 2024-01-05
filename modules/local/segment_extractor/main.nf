process SEGMENT_EXTRACTOR {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${meta.id}"
    label "few_cpu"
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'ignore' } 
    maxRetries 3
    
    input:
    tuple val(meta), path(mappedFile)

    output:
    tuple val(meta), path("segments/*_mapped-segment.tsv"), emit: mapSeg

    script:
    def outFileName = "${meta.id}_mapped-segment.tsv"

    """
    mkdir segments
    bedtools merge -i ${mappedFile} -c 1,5,5,5,5,1 -delim ";" -o count,min,max,mean,median,collapse > segments/${outFileName}
    """
}