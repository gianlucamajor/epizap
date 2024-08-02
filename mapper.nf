

include { IDX_BUILDER } from './modules/local/bowtie2/build_idx/main.nf'
include { MAPPER } from './modules/local/bowtie2/mapping/main.nf'
include { PRIMARY_MAPPED_EXTRACTOR } from './modules/local/samtools/primaryMapping/main.nf'

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
    PRIMARY_MAPPED_EXTRACTOR(MAPPER.out.mapped) 

    emit: 
    PRIMARY_MAPPED_EXTRACTOR.out

}