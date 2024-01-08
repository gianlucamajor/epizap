
process PRIMARY_MAPPED_EXTRACTOR {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    label 'few_cpu'
    tag "${meta.id}"

    input:
    tuple val(meta), path(mappingFile)
    
    output:
    path("only-primary-mapping/*_opm.{bam,sam}"), emit: opm

    script: 
    """
    echo "${meta.id} ${mappingFile}"
    mkdir only-primary-mapping
    samtools view -b -F 256 -F 4 ${mappingFile} > only-primary-mapping/${meta.id}"_opm.bam"
    """
}
