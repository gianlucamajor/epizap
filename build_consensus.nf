
params.mappingFiles = "${projectDir}/results/only-primary-mapping/chagasic_patients_opm*.bam"
params.peptides = "${projectDir}/dataSet/peptides/kDNA-aa-in-frame.fasta"
params.outdir = 'results'

log.info """\
    ===================================
    mapping Files: ${params.mappingFiles}
    peptiedes: ${params.peptides}
    outdir: ${params.outdir}
    ===================================
    """.stripIndent()


include { SEQUENCE_EXTRACTOR } from './modules/local/sequence_extractor/main.nf'
include { SEGMENT_EXTRACTOR } from './modules/local/segment_extractor/main.nf'
include { MSA } from './modules/local/msa/main.nf'
include { HMM_PROFILE_BUILDER } from './modules/local/hmm_profile_builder/main.nf'
include { HMM_CONSENSUS_BUILDER } from './modules/local/hmm_consensus_builder/main.nf'
include { LONELY } from './modules/local/lonely/main.nf'

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


    segmentsCh = SEGMENT_EXTRACTOR(mappingFilesCh)
    pepSegmentsCh = SEQUENCE_EXTRACTOR(segmentsCh, peptidesCh)


    // // mappingFilesCh.first().view()
    pepSegmentsCh.flatten()
        .filter({ it.countFasta() > 1 && it.countFasta() < 30000})
        .map{it ->
            meta = [id: it.name.replaceFirst(".fasta", ""), records: it.countFasta()]
            [meta, it]
        }
        .set{pepSegmentsChFlt}

    // pepSegmentsChFlt.view()

    pepSegmentsCh.flatten()
        .filter({ it.countFasta() < 2})
        .map{it ->
            meta = [id: it.name.replaceFirst(".fasta", ""), records: it.countFasta()]
            [meta, it]
        }
        .set{lonelyPepSegmentsChFlt}

    
    MSA(pepSegmentsChFlt) | HMM_PROFILE_BUILDER | HMM_CONSENSUS_BUILDER
    LONELY(lonelyPepSegmentsChFlt)

}
