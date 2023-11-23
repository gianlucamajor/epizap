
params.mapingFiles = "${projectDir}/results/*.bam"
params.outdir = 'results'

log.info """\
    ===================================
    mapping files: ${params.mapingFiles}
    outdir: ${params.outdir}
    ===================================
    """.stripIndent()

process primaryMappedExtractor {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${meta.id}"

    input:
    tuple val(meta), path(mappingFile)

    output:
    output:
    path("only-primary-mapping/*_opm.{bam,sam}"), emit: opm

    script: 
    """
    echo "${meta.id} ${mappingFile}"
    mkdir only-primary-mapping
    samtools view -b -F 256 -F 4 ${mappingFile} > only-primary-mapping/${meta.id}"_opm.bam"
    """
}

workflow{
    Channel.fromPath(params.mapingFiles)
        .map{ it -> 
            meta = [id:it.simpleName] 
            [meta, it]}
        .set{mappingFilesCh}
    
    primaryMappedExtractor(mappingFilesCh)
}
