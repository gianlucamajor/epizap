params.pepSeg = "${projectDir}/results/peptides-segment/*.fasta"
params.outdir = 'results'

log.info """\
    ===================================
    Peptiede Segment: ${params.pepSeg}
    outdir: ${params.outdir}
    ===================================
    """.stripIndent()

include { LONELY } from './modules/local/lonely/main.nf'

workflow{

    Channel.fromPath(params.pepSeg)
        .filter({ it.countFasta() < 2})
        .collect()
        .set{lonelyPepSegCh}
    
    LONELY(lonelyPepSegCh) | view()

}
