process SEGMENT_EXTRACTOR {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${meta.id}"
    label "med_cpu_high_memory"
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}
    maxRetries 3
    
    input:
    tuple val(meta), path(mappedFile)

    output:
    tuple val(meta), path("segments/*_mapped-segment.tsv"), emit: mapSeg

    script:
    def outFileName = "${meta.id}_mapped-segment.tsv"

    """
    mkdir segments
    bedtools merge -i ${mappedFile} -s -c 1,5,5,5,5,1 -delim ";" -o count,min,max,mean,median,collapse > segments/${outFileName}
    """
}

/*
* This process uses bedtools to intersect two files: `segmentFile` and `annotationFile`.
*       -wb	Write the original entry in B for each overlap. Useful for knowing what A overlaps.
*       for more info about bedtools intersect: https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
* The intersected output is then piped to the `cut` command to select specific columns.
*       Especially, to avoid to duplicate columns with reads info. 
* The selected columns are saved to a file in the `segments-annotated` directory with the name format `${meta.id}-annotated.tsv`.
*
* Parameters:
* - segmentFile: The file containing the segments with peptides (epitopes) mapped.
* - annotationFile: The file containing the annotations (i.g: proteins) to intersect with the segments.
* - meta.id: A identifier used to name the output file.
*
* Output:
* - A TSV file containing the segments with  intersected and selected columns from the input files.
*/
process SEGMENT_NOTE_TAKER {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${meta.id}"
    label "few_cpu"
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    tuple val(meta), path(segmentFile)
    path(annotationFile)

    output:
    tuple val(meta), path("segments-annotated/*annotated.tsv"), emit: segAnnotated
    script:
    """
    mkdir segments-annotated
    bedtools intersect -a ${segmentFile}  -b ${annotationFile} -wb | cut -f1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19 > segments-annotated/${meta.id}-annotated.tsv
    """
 
    stub:
    """
    mkdir segments-annotated
    touch segments-annotated/${meta.id}stub_annotated.tsv
    """
}