params.graph = ""
params.outdir = "results/ept-annot-reported"
params.single_reads = false
params.proteome_hits = ""
params.iedb = false
params.iedb_tcruzi_epitopes_hits = ""
params.iedb_tcruzi_epitopes = ""
params.iedb_human_epitopes_hits = ""
params.iedb_human_epitopes = ""
params.inserts_group = ""

include { EPITOPE_ANNOT_REPORTER } from './modules/local/annot_reporter/main.nf'

workflow {

    if (!params.graph) {
        error "Missing required param --graph <epitopes_graph.pickle>"
    }
    if (params.iedb && (!params.iedb_tcruzi_epitopes_hits || !params.iedb_tcruzi_epitopes)) {
        error "--iedb requires --iedb_tcruzi_epitopes_hits and --iedb_tcruzi_epitopes"
    }

    log.info """\
        ===================================
        Epizap Annotation Report
        ===================================
        Graph: ${params.graph}
        Outdir: ${params.outdir}
        Single Reads: ${params.single_reads}
        Proteome Hits: ${params.proteome_hits}
        IEDB: ${params.iedb}
        IEDB T. cruzi Epitopes Hits: ${params.iedb_tcruzi_epitopes_hits}
        IEDB T. cruzi Epitopes: ${params.iedb_tcruzi_epitopes}
        IEDB Human Epitopes Hits: ${params.iedb_human_epitopes_hits}
        IEDB Human Epitopes: ${params.iedb_human_epitopes}
        Inserts Group: ${params.inserts_group}
        ===================================
        """.stripIndent()

    Channel.fromPath(params.graph, checkIfExists: true)
        .map{ it -> [[id: it.simpleName], it] }
        .set{ graphCh }

    Channel.fromPath(params.inserts_group ?: "NO_FILE").first().set{ insertsGroupCh }
    Channel.fromPath(params.proteome_hits ?: "NO_FILE").first().set{ proteomeHitsCh }
    Channel.fromPath(params.iedb_tcruzi_epitopes_hits ?: "NO_FILE").first().set{ iedbTcruziEpitopesHitsCh }
    Channel.fromPath(params.iedb_tcruzi_epitopes ?: "NO_FILE").first().set{ iedbTcruziEpitopesCh }
    Channel.fromPath(params.iedb_human_epitopes_hits ?: "NO_FILE").first().set{ iedbHumanEpitopesHitsCh }
    Channel.fromPath(params.iedb_human_epitopes ?: "NO_FILE").first().set{ iedbHumanEpitopesCh }

    EPITOPE_ANNOT_REPORTER(
        graphCh,
        insertsGroupCh,
        proteomeHitsCh,
        params.iedb,
        iedbTcruziEpitopesHitsCh,
        iedbTcruziEpitopesCh,
        iedbHumanEpitopesHitsCh,
        iedbHumanEpitopesCh,
        params.single_reads
    )
}
