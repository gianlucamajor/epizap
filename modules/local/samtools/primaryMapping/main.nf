
process PRIMARY_MAPPED_EXTRACTOR {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    label 'many_cpu'
    tag "${meta.id}"

    input:
    tuple val(meta), path(mappingFile)
    
    output:
    tuple val(meta), path("only-primary-mapping/*_opm.{bam,sam}"), emit: opm

    script: 
    """
    echo "${meta.id} ${mappingFile}"
    mkdir only-primary-mapping
    samtools view --threads $task.cpus -b -F 256 -F 4 ${mappingFile} > only-primary-mapping/${meta.id}"_opm.bam"
    """
}
