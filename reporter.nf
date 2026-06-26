include {EPITOPE_REPORTER} from './modules/local/reporter/main.nf'

workflow epitopeReporter {
    take:
    graphCh
    insertsGroupFp

    main:
    graphCh
        .map{ it ->
            def metaGraph = [id:it.simpleName]
            [metaGraph, it]
        }.set{graphForReportCh}

    if (insertsGroupFp) {
        channel.fromPath(insertsGroupFp).first().set{insertsGroupCh}
    } else {
        channel.fromPath("NO_FILE").set{insertsGroupCh}
    }

    EPITOPE_REPORTER(graphForReportCh, insertsGroupCh)

    emit:
    jsonReport = EPITOPE_REPORTER.out.jsonReport
    fastaReport = EPITOPE_REPORTER.out.fastaReport
    epitopeGraph = EPITOPE_REPORTER.out.graphReport
    
}