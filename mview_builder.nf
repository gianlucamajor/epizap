params.baseDir = "/home/gianluca/workspace/epizap/results_04_10_2024"
params.msaFiles = "${params.baseDir}/msa/*.msc"
params.outdir = "${params.baseDir}"

log.info """\
========================================================
MSA Files = ${params.msaFiles}
========================================================
""".stripIndent()

process mview {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${meta.id}"

    input:
    tuple val(meta), path(msaFile)

    output:
    tuple val(meta), path("mview/*.html")

    script:
    def outFileName = "${meta.id}.html"
    """
    mkdir mview
    mview -in fasta -html head -css on -coloring consensus -threshold 100 -consensus on ${msaFile} > mview/${outFileName}
    """
}
workflow{
    Channel.fromPath(params.msaFiles)
    .map{ it ->
        fnTokens = it.fileName.toString() - "control_and_chagasic_patients_" //remove prefix string
        fnTokens = fnTokens - "-msa.msc" // remove sufix string
        fnTokens = fnTokens.replace("_", "-")
        metaMap = [id:fnTokens]
        [metaMap, it]
    }.set{msas}

    mview(msas)
}