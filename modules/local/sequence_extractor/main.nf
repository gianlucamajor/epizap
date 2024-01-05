process SEQUENCE_EXTRACTOR {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${metaSeg.id}"
    label "few_cpu"
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'ignore' } 
    maxRetries 3

    input:
    tuple val(metaSeg), path(segmentFile)
    tuple val(metaPep), path(peptideFile)

    output:
    path("peptides-segment/*.fasta")

    script:
    """
    mkdir peptides-segment
    peptide_extractor.py  ${segmentFile} ${peptideFile} -o peptides-segment -pf ${metaSeg.id}
    """
}
