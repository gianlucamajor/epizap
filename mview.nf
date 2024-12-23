params.baseDir = "/home/gianluca/workspace/epizap/results"
params.msaFiles = "${params.baseDir}/msa/*.msc"
params.outdir = "${params.baseDir}"


include { MVIEW } from './modules/local/mview/main.nf'

workflow mview {
    take:
    msaChn

    main:
    msaChn.map{ it ->
        def fnTokens = it.fileName.toString() - "-msa.msc"
        def metaMap = [id:fnTokens]
        [metaMap, it]
    }.set{msas}

    MVIEW(msas)

    emit:
    MVIEW.out
}