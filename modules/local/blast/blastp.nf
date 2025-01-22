params.sequences = "${projectDir}/dataSet/epitopes/epitopes-firt-10-mmc2.faa"
params.db = "${projectDir}/blast/db/"
params.dbName = "epitope_table_04_08_22_T_Cruzi_1659634014"
params.outDir = 'blast'

process blastp {
    publishDir "${params.outDir}", mode: 'copy', overwrite: true
    tag "${meta.id}"
    label "many_cpu"

    input:
    tuple val(meta), path(sequences)
    tuple val(metaDb), path(db)

    output:
    tuple val(meta), path("${sequences.simpleName}*.tsv"), emit: blast_results

    
    script:
    """
    blastp -outfmt '6 std ssciname sseqid sgi' \\
    -num_threads ${task.cpus} \\
    -db ${db}/${metaDb.id} \\
    -query ${sequences} > ${sequences.simpleName}.tsv
    """
}

workflow {

    log.info """\
    ===================================
    sequences: ${params.sequences}
    database name: ${params.dbName}
    database: ${params.db}
    outdir: ${params.outDir}
    ===================================
    """.stripIndent()
    
    Channel.fromPath(params.db)
        .map{it ->
            def meta = [id: params.dbName]
            [meta, it]
        }
        .first()
        .set{dbCh}

    Channel.fromPath(params.sequences)
        .map{it ->
            def meta = [id: it.simpleName]
            [meta, it]
        }
        .set{sequencesCh}

    blastp(sequencesCh, dbCh)
}