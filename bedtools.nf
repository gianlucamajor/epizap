params.mappingFiles = "${projectDir}/results/only-primary-mapping/*.bam"
params.outdir = 'results'

log.info """\
    ===================================
    reference: ${params.mappingFiles}
    outdir: ${params.outdir}
    ===================================
    """.stripIndent()

process segmentExtractor {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${meta.id}"

    input:
    tuple val(meta), path(mappedFile)

    output:
    tuple val(meta), path("*_mapped-segment.tsv"), emit: mapSeg

    script:
    def outFileName = "${meta.id}_mapped-segment.tsv"

    """
    bedtools merge -i ${mappedFile} -c 1,5,5,5,5,1 -delim ";" -o count,min,max,mean,median,collapse > ${outFileName} 
    """

}

process peptide_extractor {
    tag "${meta.id}"
    
    input:
    tuple val(meta), path(segmentFile)

    output:
    stdout

    script:
    """
    python3 ${projectDir}/bin/peptide_extractor.py  ${segmentFile}
    """

}

workflow{

    Channel.fromPath(params.mappingFiles)
        .map{it -> 
            meta = [id: it.simpleName]
            [meta, it]
        }       
        .set{mappingFilesCh}

    segmentExtractor(mappingFilesCh) | peptide_extractor

}
