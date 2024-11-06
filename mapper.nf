

include { IDX_BUILDER } from './modules/local/bowtie2/build_idx/main.nf'
include { MAPPER } from './modules/local/bowtie2/mapping/main.nf'
include { MAX_SCORE_MAPPING_EXTRACTOR } from './modules/local/samtools/max_score_mapping/main.nf'
// include { PRIMARY_MAPPED_EXTRACTOR } from './modules/local/samtools/primaryMapping/main.nf' # TO BE REMOVED

workflow mapper {
    take:
    ref
    reads
    
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

    IDX_BUILDER(inputRefIdx)
    MAPPER(IDX_BUILDER.out.idx, readsCh)
    MAX_SCORE_MAPPING_EXTRACTOR(MAPPER.out.mapped)
    // PRIMARY_MAPPED_EXTRACTOR(MAPPER.out.mapped) # TO BE REMOVED

    emit: 
    MAX_SCORE_MAPPING_EXTRACTOR.out

}