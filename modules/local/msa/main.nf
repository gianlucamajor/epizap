process MSA {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "${pepSegFile.name}"
    label "med_cpu"
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'ignore' } 
    maxRetries 3

    input:
    tuple val(pepSegMeta), path(pepSegFile)

    output:
    path("msa/*-msa.msc")

    script:
    """
    mkdir msa
    muscle -align ${pepSegFile} -output msa/${pepSegFile.baseName}-msa.msc -threads ${task.cpus}
    """
}