params.reads = "${projectDir}/dataSet/k-dna-insert-reads/*.fastq"
params.peptides = "${projectDir}/dataSet/peptides/kDNA-aa-in-frame.fasta"
params.ref = "${projectDir}/dataSet/ref/y_strain_minicircles.fasta"
params.outdir = "results"
params.mapper = ""
params.annotation = ""
params.graph_segment_threshold = 1.00


params.pepSeg = "${projectDir}/${params.outdir}/peptides-segment/*.fasta"

include { mapper } from "./mapper"
include { segmentRetriever } from "./segment_retriever"
include { segmentGraphComputer } from './segment_graph_computer'
include { msa } from "./msa"
include { mview } from "./mview"
include { predicator } from "./predicator"
include { lonelyPeptideRetriever } from './lonely_peptide_retriever'
include { msaCoreAndLonelyEpitopesJoiner } from './msa_core_and_lonely_epitopes_joiner'

workflow {

log.info """\
    ===================================
    Reads: ${params.reads}
    Peptides: ${params.peptides}
    Reference: ${params.ref}
    Annotation: ${params.annotation}
    Graph Segment threshold: ${params.graph_segment_threshold}
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

    segmentRetriever(mapping, params.peptides, params.annotation)
        .set{segmentRetriverResult}
    msa(segmentRetriverResult.peptidesFromSegment)
    segmentGraphComputer(segmentRetriverResult.segments)

    mview(msa.out)
    predicator(mview.out)


    lonelyPeptideRetriever(segmentRetriverResult.peptidesFromSegment)

    lonelyPeptideRetriever.out.mix(predicator.out)
        .set{epitopesByMSACoreAndLonely}
        
    msaCoreAndLonelyEpitopesJoiner(epitopesByMSACoreAndLonely)
    
}