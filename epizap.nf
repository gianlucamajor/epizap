params.reads = "${projectDir}/dataSet/k-dna-insert-reads/*.fastq"
params.peptides = "${projectDir}/dataSet/peptides/kDNA-aa-in-frame.fasta"
params.ref = "${projectDir}/dataSet/ref/y_strain_minicircles.fasta"
params.outdir = "results"
params.mapper = ""
params.annotation = "${projectDir}/dataSet/ref/annotation_genomic_dummy.gff" //.gff
params.graph_segment_threshold = 1.00
params.inserts_group = ""


params.pepSeg = "${projectDir}/${params.outdir}/peptides-segment/*.fasta"

include { mapper } from "./mapper"
include { segmentRetriever } from "./segment_retriever"
include { segmentGraphComputer } from './segment_graph_computer'
include { peptideClusteringByCC } from './peptide_clustering_by_CC'
include { msa } from "./msa"
include { mview } from "./mview"
include { predicator } from "./predicator"
include { graphUpdater } from './graph_updater'
include { lonelyPeptideRetriever } from './lonely_peptide_retriever'
include { msaCoreAndLonelyEpitopesJoiner } from './msa_core_and_lonely_epitopes_joiner'
include { epitopeReporter } from './reporter'    

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
    Inserts Group: ${params.inserts_group}

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

    //  params.annotation shoudl be a optional file
    segmentRetriever(mapping, params.peptides, params.annotation)
        .set{segmentRetriverResult}
    segmentGraphComputer(segmentRetriverResult.segments)

    
    // call peptide clustering graph with graph as parameter.
    peptideClusteringByCC(segmentGraphComputer.out.segmentGraph, params.peptides)

    msa(peptideClusteringByCC.out.listOfPeptidesByCC)
    mview(msa.out)
    predicator(mview.out)
    
    graphUpdater(peptideClusteringByCC.out.graph, predicator.out.msaEpitopeReport, segmentRetriverResult.aFeatures)
    
    epitopeReporter(graphUpdater.out.graph, params.inserts_group)

    // lonelyPeptideRetriever(segmentRetriverResult.peptidesFromSegment)

    // lonelyPeptideRetriever.out.mix(predicator.out)
    //     .set{epitopesByMSACoreAndLonely}        
    // msaCoreAndLonelyEpitopesJoiner(epitopesByMSACoreAndLonely)
    
}