params.sequences = "${projectDir}/dataSet/iedb/epitope_table_04_08_22_T_Cruzi_1659634014.fasta"
params.dbType = 'prot'
params.outDir = 'blast/db'



process make_blast_db {
    publishDir "${params.outDir}", mode: 'copy', overwrite: true
    tag "${meta.id}"

    input:
    tuple val(meta), path(sequences)

    output:
    tuple val(meta), path("${sequences.simpleName}*"), emit: blast_db

    
    script:
    """
    makeblastdb -in ${sequences} -out ${sequences.simpleName} -dbtype ${params.dbType}
    """
}

workflow {
    log.info """\
    ===================================
    sequences: ${params.sequences}
    outdir: ${params.outDir}
    ===================================
    """.stripIndent()

    Channel.fromPath(params.sequences)
        .map{it ->
            def meta = [id: it.simpleName]
            [meta, it]
        }
        .set{sequencesCh}

    make_blast_db(sequencesCh)
}