params.outdir = "results_graph_21_01_2025/epitopes-by-msa-core-and-lonely/results-msa-snippet/mafft"
// params.msaFiles = "results_graph_21_01_2025/epitopes-by-msa-core-and-lonely/epitopes-by-cc/*.fasta"
// params.msa_in_files = "results-msa-snippet/msa/*.msc"
params.msaFiles = "results_graph_21_01_2025/epitopes-by-msa-core-and-lonely/cc-clustered/*.fasta"
params.msa_in_files = "results_graph_21_01_2025/epitopes-by-msa-core-and-lonely/results-msa-snippet/mafft/mafft-output/*.mafft"



include { MSA } from './modules/local/msa/main.nf'
include { MVIEW } from './modules/local/mview/main.nf'


workflow  {
    log.info """\
    ===================================
    outdir: ${params.outdir}
    msa files: ${params.msaFiles}
    ===================================
    """.stripIndent()
    Channel.fromPath(params.msaFiles)
            .map{ it ->
                def metaAlign = [id: it.name.replaceFirst(".fasta", ""), records: it.countFasta(), algorithm: "super5"]
                [metaAlign, it]
            }.set{msaPep}

    Channel.fromPath(params.msa_in_files)
    .map{ it ->
        def fnTokens = it.fileName.toString() - "-msa.msc"
        def metaMap = [id:fnTokens]
        [metaMap, it]
    }.set{msas}

    // MSA(msaPep)
    MVIEW(msas)    
}