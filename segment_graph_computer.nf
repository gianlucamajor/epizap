include { SEGMENT_GRAPH_COMPUTER } from './modules/local/segment_graph_computer/main.nf'

workflow segmentGraphComputer{
    take:
    segments
    
    main:    
    SEGMENT_GRAPH_COMPUTER(segments)
   
    
    emit:
    segmentGraph = SEGMENT_GRAPH_COMPUTER.out

}
