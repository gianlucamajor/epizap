process LONELY {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    tag " Number of segments with lonely peptide: ${files.size()}"
    label "med_cpu"
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'ignore' } 
    maxRetries 3

    input:
    path(files)

    output:
    path("lonely/lonely_peptides.fasta")


    script:
    """
    mkdir lonely
    lonely_peptide_generator.py ${files} -o lonely/lonely_peptides.fasta --threads ${task.cpus}
    """ 
}