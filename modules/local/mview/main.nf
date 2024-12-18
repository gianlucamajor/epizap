process MVIEW {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${meta.id}"

    input:
    tuple val(meta), path(msaFile)

    output:
    path("mview/*.html")
    
    script:
    def outFileName = "${meta.id}.html"
    """
    mkdir mview
    mview -in fasta -html head -css on -coloring consensus -threshold 100 -consensus on ${msaFile} > mview/${outFileName}
    """
}