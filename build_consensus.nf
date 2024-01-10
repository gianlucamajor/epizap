
params.pepSeg = "${projectDir}/results/peptides-segment/*.fasta"
params.outdir = 'results'

log.info """\
    ===================================
    Peptiede Segment: ${params.pepSeg}
    outdir: ${params.outdir}
    ===================================
    """.stripIndent()

include { MSA } from './modules/local/msa/main.nf'
include { HMM_PROFILE_BUILDER } from './modules/local/hmm_profile_builder/main.nf'
include { HMM_CONSENSUS_BUILDER } from './modules/local/hmm_consensus_builder/main.nf'

workflow{

    Channel.fromPath(params.pepSeg)
        .set{pepSegCh}


    pepSegCh
        .filter({ it.countFasta() > 1 && it.countFasta() < 30000})
        .map{it ->
            meta = [id: it.name.replaceFirst(".fasta", ""), records: it.countFasta()]
            [meta, it]
        }
        .set{msaPepSegCh}    
    
    MSA(msaPepSegCh) | HMM_PROFILE_BUILDER | HMM_CONSENSUS_BUILDER

}
