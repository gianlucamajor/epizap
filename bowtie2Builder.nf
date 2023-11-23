params.ref = "${projectDir}/dataSet/ref/y_strain_minicircles.fasta"
params.outdir = 'results'

log.info """\
    ===================================
    reference: ${params.ref}
    outdir: ${params.outdir}
    ===================================
    """.stripIndent()



process indexBuilder {
    publishDir "${params.outdir}", mode: 'copy', overwrite: false
    tag "${meta}"
    cpus 6
    input:
    tuple val(meta), path(fastaRef)

    output:
    path("${fastaRef.simpleName}-idx"), emit: idx

    script:
    """
    mkdir "${fastaRef.simpleName}-idx"
    bowtie2-build --threads $task.cpus ${fastaRef} "${fastaRef.simpleName}-idx/${fastaRef.simpleName}"
    """
}

workflow{
    Channel.fromPath(params.ref)
        .map{it -> [it.simpleName, it]}
        .groupTuple()
        .set{inputRefIdx}

    indexBuilder(inputRefIdx)

}
