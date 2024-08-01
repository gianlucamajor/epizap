process LONELY {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag "number of segments with lonely peptide: ${files.size()}"
    label "med_cpu"
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    path(files)

    output:
    path("lonely/lonely_peptides.fasta"), emit: lonelyPep


    script:
    """
    mkdir lonely
    lonely_peptide_generator.py ${files} -o lonely/lonely_peptides.fasta --threads ${task.cpus}
    """ 
}