params.reads = "${projectDir}/dataSet/k-dna-insert-reads/*.fastq"
params.peptides = "${projectDir}/dataSet/peptides/kDNA-aa-in-frame.fasta"
params.ref = "${projectDir}/dataSet/ref/y_strain_minicircles.fasta"
params.outdir = 'results'

log.info """\
    ===================================
    Reads: ${params.reads}
    Peptides: ${params.peptides}
    Reference: ${params.ref}
    outdir: ${params.outdir}
    
    ===================================
    """.stripIndent()

include { IDX_BUILDER } from './modules/local/bowtie2/build_idx/main.nf'
include { MAPPER } from './modules/local/bowtie2/mapping/main.nf'
include { PRIMARY_MAPPED_EXTRACTOR } from './modules/local/samtools/primaryMapping/main.nf'
include { SEGMENT_EXTRACTOR } from './modules/local/segment_extractor/main.nf'
include { SEQUENCE_EXTRACTOR } from './modules/local/sequence_extractor/main.nf'


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
    
    Channel.fromPath(params.peptides)
        .map{it ->
            meta = [id: it.simpleName]
            [meta, it]
        }
        .first()
        .set{peptidesCh}

    IDX_BUILDER(inputRefIdx)
    MAPPER(IDX_BUILDER.out.idx, readsCh)
    PRIMARY_MAPPED_EXTRACTOR(MAPPER.out.mapped) 
    SEGMENT_EXTRACTOR(PRIMARY_MAPPED_EXTRACTOR.out.opm)
    SEQUENCE_EXTRACTOR(SEGMENT_EXTRACTOR.out.mapSeg, peptidesCh)

}