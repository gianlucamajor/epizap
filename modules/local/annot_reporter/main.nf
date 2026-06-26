process EPITOPE_ANNOT_REPORTER {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${meta.id}"
    label "med_cpu"
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    tuple val(meta), path(graphFile)
    path(insertsGroupFile)
    path(proteomeHitsFile)
    val(iedb)
    path(iedbTcruziEpitopesHitsFile)
    path(iedbTcruziEpitopesFile)
    path(iedbHumanEpitopesHitsFile)
    path(iedbHumanEpitopesFile)
    val(singleReads)

    output:
    path("epitopes_report/*.json"), emit: jsonReport
    path("epitopes_report/*.fasta"), emit: fastaReport
    path("epitopes_report/*.pickle"), emit: graphReport

    script:
    def insertsGroupArg = insertsGroupFile.name != 'NO_FILE' ? "--inserts-group ${insertsGroupFile}" : ""
    def proteomeHitsArg = proteomeHitsFile.name != 'NO_FILE' ? "--proteome-hits ${proteomeHitsFile}" : ""
    def iedbArg = iedb ? "--iedb" : ""
    def iedbTcruziEpitopesHitsArg = iedbTcruziEpitopesHitsFile.name != 'NO_FILE' ? "--iedb-tcruzi-epitopes-hits ${iedbTcruziEpitopesHitsFile}" : ""
    def iedbTcruziEpitopesArg = iedbTcruziEpitopesFile.name != 'NO_FILE' ? "--iedb-tcruzi-epitopes ${iedbTcruziEpitopesFile}" : ""
    def iedbHumanEpitopesHitsArg = iedbHumanEpitopesHitsFile.name != 'NO_FILE' ? "--iedb-human-epitopes-hits ${iedbHumanEpitopesHitsFile}" : ""
    def iedbHumanEpitopesArg = iedbHumanEpitopesFile.name != 'NO_FILE' ? "--iedb-human-epitopes ${iedbHumanEpitopesFile}" : ""
    def singleReadsArg = singleReads ? "--single-reads" : ""
    """
    mkdir -p epitopes_report
    epitope_reporter.py ${graphFile} -o epitopes_report/ ${insertsGroupArg} ${proteomeHitsArg} ${iedbArg} ${iedbTcruziEpitopesHitsArg} ${iedbTcruziEpitopesArg} ${iedbHumanEpitopesHitsArg} ${iedbHumanEpitopesArg} ${singleReadsArg}
    cp ${graphFile} epitopes_report/
    """
}
