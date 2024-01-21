process MAPPER {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${metaReads.id}"
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}

    input:
    tuple val(metaParent), path(bwt2Idx)
    tuple val(metaReads), path(reads)
    label 'many_cpu'

    output:
    tuple val(metaReads), path("mapping/*.{bam,sam}"), emit: mapped

    script:
    def bwt2IdxPath = "${bwt2Idx}/${metaParent.id}"

    """
    mkdir mapping
    bowtie2 -a --threads $task.cpus -x ${bwt2IdxPath} -U ${reads} | samtools view -bS - | samtools sort -o mapping/${metaReads.id}.bam
    """

}
