

include { IDX_BUILDER } from './modules/local/bowtie2/build_idx/main.nf'
include { MAPPER } from './modules/local/bowtie2/mapping/main.nf'
include { PRIMARY_MAPPED_EXTRACTOR } from './modules/local/samtools/primaryMapping/main.nf'
include { SEGMENT_EXTRACTOR } from './modules/local/segment_extractor/main.nf'
include { SEQUENCE_EXTRACTOR } from './modules/local/sequence_extractor/main.nf'


workflow mapper {
    take:
    ref
    reads
    peptides

    
    main:
    Channel.fromPath(ref)
        .map{ it -> 
                metaRef = [id:it.simpleName]
                [metaRef, it]
        }
        .first()
        .set{inputRefIdx}

        Channel.fromPath(reads)
        .map{ it -> 
            metaReads = [id:it.simpleName]
            [metaReads, it]
        }
        .set{readsCh}
    
    Channel.fromPath(peptides)
        .map{it ->
            metaPep = [id:it.simpleName]
            [metaPep, it]
        }
        .first()
        .set{peptidesCh}


    IDX_BUILDER(inputRefIdx)
    MAPPER(IDX_BUILDER.out.idx, readsCh)
    PRIMARY_MAPPED_EXTRACTOR(MAPPER.out.mapped) 
    SEGMENT_EXTRACTOR(PRIMARY_MAPPED_EXTRACTOR.out.opm)
    SEQUENCE_EXTRACTOR(SEGMENT_EXTRACTOR.out.mapSeg, peptidesCh)

    emit: 
    SEQUENCE_EXTRACTOR.out

}