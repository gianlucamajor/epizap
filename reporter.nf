include {EPITOPE_REPORTER} from './modules/local/reporter/main.nf'

workflow epitopeReporter {
    take:
    graphCh

    main:
    graphCh
        .map{ it -> 
            def metaGraph = [id:it.simpleName]
            [metaGraph, it]
        }.set{graphForReportCh}
    
    EPITOPE_REPORTER(graphForReportCh)

    emit:
    jsonReport = EPITOPE_REPORTER.out.jsonReport
    fastaReport = EPITOPE_REPORTER.out.fastaReport
    epitopeGraph = EPITOPE_REPORTER.out.graphReport
    
}