params.bwt2Idx = "${projectDir}/results/y_strain_minicircles-idx/*"
params.reads = "${projectDir}/dataSet/k-dna-insert-reads/*.fastq"
params.outdir = 'results'

log.info """\
    ===================================
    bwt2Idx: ${params.bwt2Idx}
    outdir: ${params.outdir}
    reads: ${params.reads}
    ===================================
    """.stripIndent()

process mapper {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${metaReads}"

    input:
    tuple val(metaParent), val(metaIdx), path(bwt2Idx)
    tuple val(metaReads), path(reads) 
    
    cpus 6
    
    output:
    tuple val(metaReads), path("*.{bam,sam}"), emit: mapped
    
    script:
    def bwt2IdxPath = "${metaParent}/${metaIdx[0]}"

    """
    bowtie2 -a --threads $task.cpus -x ${bwt2IdxPath} -U ${reads} | samtools view -bS - | samtools sort -o ${metaReads}.bam
    """
}


workflow{

    Channel.fromPath(params.bwt2Idx)
        .map(it -> [it.getParent(), it.simpleName, it])
        .groupTuple()
        .first()
        .set{bwt2IdxCh}

    Channel.fromPath(params.reads)
        .map{it -> [it.simpleName, it]}
        .groupTuple()
        .set{readsCh}


    mapper(bwt2IdxCh, readsCh)

}
