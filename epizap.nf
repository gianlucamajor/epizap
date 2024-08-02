params.reads = "${projectDir}/dataSet/k-dna-insert-reads/*.fastq"
params.peptides = "${projectDir}/dataSet/peptides/kDNA-aa-in-frame.fasta"
params.ref = "${projectDir}/dataSet/ref/y_strain_minicircles.fasta"
params.outdir = "results"
params.mapper = ""

params.pepSeg = "${projectDir}/${params.outdir}/peptides-segment/*.fasta"

log.info """\
    ===================================
    Reads: ${params.reads}
    Peptides: ${params.peptides}
    Reference: ${params.ref}
    Outdir: ${params.outdir}
    Peptiede Segment: ${params.pepSeg}
    Mapper file: ${params.mapper}
    
    ===================================
    """.stripIndent()


include { mapper } from "./mapper"
include { segmentRetriever } from "./segment_retriever"
include { consensusBuilder } from './consensus_builder'
include { lonelyPeptideRetriever } from './lonely_peptide_retriever'
include { preEpitopeRetriever } from './pre_epitope_retriever'

workflow {

    if( ! params.mapper){
        mapper(params.ref, params.reads)
        mapping = mapper.out
    } else {
        Channel.fromPath(params.mapper)
            .map{ it -> 
                metaMap = [id:it.simpleName]
                [metaMap, it]
            }.set{mapping}
    }

    segmentRetriever(mapping, params.peptides)
        .set{peptidesBySegments}
    
    consensusBuilder(peptidesBySegments)
    lonelyPeptideRetriever(peptidesBySegments)

    lonelyPeptideRetriever.out.mix(consensusBuilder.out)
        .set{preEpitopes}
        
    preEpitopeRetriever(preEpitopes)
    
}