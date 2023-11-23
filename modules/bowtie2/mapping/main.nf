process MAPPER {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${metaReads.id}"

    input:
    tuple val(metaParent), path(bwt2Idx)
    tuple val(metaReads), path(reads) 
    
    cpus 6
    
    output:
    tuple val(metaReads), path("*.{bam,sam}"), emit: mapped
    
    script:
    def bwt2IdxPath = "${bwt2Idx}/${metaParent.id}"

    """
    bowtie2 -a --threads $task.cpus -x ${bwt2IdxPath} -U ${reads} | samtools view -bS - | samtools sort -o ${metaReads.id}.bam
    """
    
}