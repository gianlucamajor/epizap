params.baseDir = "/home/gianluca/workspace/epizap/results"
params.msaFiles = "${params.baseDir}/msa/*.msc"
params.outdir = "${params.baseDir}"


include { MVIEW } from './modules/local/mview/main.nf'

workflow mview {
    take:
    msaChn

    main:
    msaChn.map{ it ->
        def fnTokens = it.fileName.toString() - "control_and_chagasic_patients_" //remove prefix string
        fnTokens = fnTokens - "-msa.msc" // remove sufix string
        fnTokens = fnTokens.replace("_", "-")
        def metaMap = [id:fnTokens]
        [metaMap, it]
    }.set{msas}

    MVIEW(msas)

    emit:
    MVIEW.out
}