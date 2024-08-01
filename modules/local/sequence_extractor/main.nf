process SEQUENCE_EXTRACTOR {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${metaSeg.id}"
    label "med_cpu_high_memory"
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    tuple val(metaSeg), path(segmentFile)
    tuple val(metaPep), path(peptideFile)

    output:
    path("peptides-segment/*.fasta"), emit: pepSeg

    script:
    """
    mkdir peptides-segment
    peptide_extractor.py  ${segmentFile} ${peptideFile} -o peptides-segment -pf ${metaSeg.id} -t ${task.cpus}
    """
}
