params.reads = "${projectDir}/dataSet/k-dna-insert-reads/*.fastq"
params.peptides = "${projectDir}/dataSet/peptides/kDNA-aa-in-frame.fasta"
params.ref = "${projectDir}/dataSet/ref/y_strain_minicircles.fasta"
params.outdir = "results"
params.mapper = ""

params.pepSeg = "${projectDir}/${params.outdir}/peptides-segment/*.fasta"

include { mapper } from "./mapper"
include { segmentRetriever } from "./segment_retriever"
include { msa } from "./msa"
include { mview } from "./mview"
include { predicator } from "./predicator"
include { lonelyPeptideRetriever } from './lonely_peptide_retriever'
include { msaCoreAndLonelyEpitopesJoiner } from './msa_core_and_lonely_epitopes_joiner.nf'

workflow {

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

    if( ! params.mapper){
        mapper(params.ref, params.reads)
        mapping = mapper.out
    } else {
        Channel.fromPath(params.mapper)
            .map{ it -> 
                def metaMap = [id:it.simpleName]
                [metaMap, it]
            }.set{mapping}
    }

    segmentRetriever(mapping, params.peptides)
        .set{peptidesBySegments}
    msa(peptidesBySegments)

    mview(msa.out)
    predicator(mview.out)


    lonelyPeptideRetriever(peptidesBySegments)

    lonelyPeptideRetriever.out.mix(predicator.out)
        .set{epitopesByMSACoreAndLonely}
        
    msaCoreAndLonelyEpitopesJoiner(epitopesByMSACoreAndLonely)
    
}