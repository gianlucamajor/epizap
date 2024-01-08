params.ref = "${projectDir}/dataSet/ref/y_strain_minicircles.fasta"
params.reads = "${projectDir}/dataSet/k-dna-insert-reads/*.fastq"
params.outdir = 'results'

log.info """\
    ===================================
    reference: ${params.ref}
    outdir: ${params.outdir}
    reads: ${params.reads}
    ===================================
    """.stripIndent()

include { IDX_BUILDER } from './modules/local/bowtie2/build/main.nf'
include { MAPPER } from './modules/local/bowtie2/mapping/main.nf'
include { PRIMARY_MAPPED_EXTRACTOR } from './modules/local/samtools/primaryMapping/main.nf'

workflow {

    Channel.fromPath(params.ref)
        .map{ it -> 
                meta = [id:it.simpleName]
                [meta, it]
        }
        .first()
        .set{inputRefIdx}

        Channel.fromPath(params.reads)
        .map{ it -> 
            meta = [id:it.simpleName]
            [meta, it]
        }
        .set{readsCh}

    IDX_BUILDER(inputRefIdx)
    MAPPER(IDX_BUILDER.out.idx, readsCh)
    PRIMARY_MAPPED_EXTRACTOR(MAPPER.out.mapped) 
    // search by intervals 
    //    
    
}