include { PEPTIDE_CLUSTERING_CC } from './modules/local/peptide_clustering_cc/main.nf'


workflow peptideClusteringByCC {
    take:
    graphFile
    pepFile
    
    main:
    Channel.fromPath(pepFile)
        .map{it ->
            def metaPep = [id:it.simpleName]
            [metaPep, it]
        }
        .first()
        .set{peptidesCh}

    graphFile
            .map{ it -> 
                def metaGraph = [id:it.simpleName]
                [metaGraph, it]
            }.set{graphCh}

    PEPTIDE_CLUSTERING_CC(graphCh, peptidesCh)
   
    
    emit:
    listOfPeptidesByCC = PEPTIDE_CLUSTERING_CC.out.peptidesByCC
    graph = PEPTIDE_CLUSTERING_CC.out.segmentGraph

}
