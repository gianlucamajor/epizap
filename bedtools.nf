params.mappingFiles = "${projectDir}/results/only-primary-mapping/*.bam"
params.peptides = "${projectDir}/dataSet/peptides/kDNA-aa-in-frame.fasta"
params.outdir = 'results'

log.info """\
    ===================================
    mapping Files: ${params.mappingFiles}
    peptiedes: ${params.peptides}
    outdir: ${params.outdir}
    ===================================
    """.stripIndent()

process segmentExtractor {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${meta.id}"
    label "few_cpu"

    input:
    tuple val(meta), path(mappedFile)

    output:
    tuple val(meta), path("segments/*_mapped-segment.tsv"), emit: mapSeg

    script:
    def outFileName = "${meta.id}_mapped-segment.tsv"

    """
    mkdir segments
    bedtools merge -i ${mappedFile} -c 1,5,5,5,5,1 -delim ";" -o count,min,max,mean,median,collapse > segments/${outFileName}
    """

}

process peptide_extractor {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${metaSeg.id}"
    label "med_cpu"

    input:
    tuple val(metaSeg), path(segmentFile)
    tuple val(metaPep), path(peptideFile)

    output:
    path("peptides-segment/*.fasta")

    script:
    """
    mkdir peptides-segment
    python3 ${projectDir}/bin/peptide_extractor.py  ${segmentFile} ${peptideFile} -o peptides-segment -pf ${metaSeg.id}
    """
}

process msa {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${pepSegFile.name}"
    label "many_cpu"
    memory { task.memory * task.attempt } 
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'ignore' } 
    maxRetries 3

    input:
    tuple val(pepSegMeta), path(pepSegFile)

    output:
    path("msa/*-msa.msc")

    script:
    """
    mkdir msa
    muscle -align ${pepSegFile} -output msa/${pepSegFile.baseName}-msa.msc -threads ${task.cpus}
    """
}

process hmm_builder{
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    label "few_cpu"

    input:
    path(msaFile)

    output:
    path("hmm/*.hmm")


    script:
    """
    mkdir hmm
    hmmbuild hmm/${msaFile.baseName}.hmm ${msaFile}
    """
}

process hmm_emit{
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    label "few_cpu"

    input:
    path(hmmFile)

    output:
    path("consensus/*.consensus")


    script:
    """
    mkdir consensus
    hmmemit -c ${hmmFile} > consensus/${hmmFile.baseName}.consensus
    """

}

workflow{

    Channel.fromPath(params.mappingFiles)
        .map{it ->
            meta = [id: it.simpleName]
            [meta, it]
        }
        .set{mappingFilesCh}

    Channel.fromPath(params.peptides)
        .map{it ->
            meta = [id: it.simpleName]
            [meta, it]
        }
        .first()
        .set{peptidesCh}


    segmentsCh = segmentExtractor(mappingFilesCh)
    pepSegmentsCh = peptide_extractor(segmentsCh, peptidesCh)


    // // mappingFilesCh.first().view()
    pepSegmentsCh.flatten()
        .filter({ it.countFasta() > 1 && it.countFasta() < 30000})
        .map{it ->
            meta = [id: it.name.replaceFirst(".fasta", ""), records: it.countFasta(), memory: "high_memory"]
            [meta, it]
        }
        .set{pepSegmentsChFlt}

    // pepSegmentsChFlt.view()
    
    msa(pepSegmentsChFlt) | hmm_builder | hmm_emit

}
