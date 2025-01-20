include { SEGMENT_EXTRACTOR } from './modules/local/segment_extractor/main.nf'
include { SEQUENCE_EXTRACTOR } from './modules/local/sequence_extractor/main.nf'

workflow segmentRetriever{
    take:
    mapping
    peptides
    
    main:
    Channel.fromPath(peptides)
        .map{it ->
            def metaPep = [id:it.simpleName]
            [metaPep, it]
        }
        .first()
        .set{peptidesCh}
    
    SEGMENT_EXTRACTOR(mapping)
    SEQUENCE_EXTRACTOR(SEGMENT_EXTRACTOR.out, peptidesCh)
    
    emit:
    segments = SEGMENT_EXTRACTOR.out
    peptidesFromSegment = SEQUENCE_EXTRACTOR.out
    
}
