include { GRAPH_UPDATER } from './modules/local/graph_updater/main.nf'

workflow graphUpdater {
    take:
    graph
    report
    aFeaturesCh

    main:

    graph
    .map{ it -> 
            def metaGraph = [id:it.simpleName]
            [metaGraph, it]
        }.set{graphCh}
    
    report
        .map{ it -> 
            def metaReport = [id:it.simpleName]
            [metaReport, it]
        }.set{reportCh}



    GRAPH_UPDATER(graphCh, reportCh, aFeaturesCh)
    aFeaturesCh.view()

    
    emit:
    graph = GRAPH_UPDATER.out.graph
    
}