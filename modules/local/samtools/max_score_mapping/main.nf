
process MAX_SCORE_MAPPING_EXTRACTOR {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    label 'one_cpu'
    tag "${meta.id}"
    errorStrategy {task.attempt <= 3 ? 'retry' : 'ignore'}

    input:
    tuple val(meta), path(mapping_file)
    
    output:
    tuple val(meta), path("max-score-mapping/*_msm_sorted.{bam,sam}"), emit: msm 

    script: 
    """
    echo "${meta.id} ${mapping_file}"
    mkdir max-score-mapping
    bam_handler.py ${mapping_file} max-score-mapping
    samtools sort max-score-mapping/${meta.id}_msm.bam -o max-score-mapping/${meta.id}_msm_sorted.bam
    """
}
