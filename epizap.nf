params.reads = "${projectDir}/dataSet/k-dna-insert-reads/*.fastq"
params.peptides = "${projectDir}/dataSet/peptides/kDNA-aa-in-frame.fasta"
params.ref = "${projectDir}/dataSet/ref/y_strain_minicircles.fasta"
params.outdir = "results"

params.pepSeg = "${projectDir}/${params.outdir}/peptides-segment/*.fasta"

log.info """\
    ===================================
    Reads: ${params.reads}
    Peptides: ${params.peptides}
    Reference: ${params.ref}
    Outdir: ${params.outdir}
    Peptiede Segment: ${params.pepSeg}
    
    ===================================
    """.stripIndent()


include { mapper } from "./mapper"
include { consensusBuilder } from './consensus_builder'
include { lonelyPeptideRetriever } from './lonely_peptide_retriever'
include { preEpitopeRetriever } from './pre_epitope_retriever'

workflow {
    mapper(params.ref, params.reads, params.peptides)
    consensusBuilder(mapper.out)
    lonelyPeptideRetriever(mapper.out)

    lonelyPeptideRetriever.out.mix(consensusBuilder.out)
        .set{preEpitopes}
        
    preEpitopeRetriever(preEpitopes)
    
}