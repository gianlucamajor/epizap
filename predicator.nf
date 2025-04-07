include { PREDICATOR } from './modules/local/predicator/main.nf'

workflow predicator{ 
    take:
    mviewCh
    
    main:
    mviewCh
        .flatten()
        .collect()
        .set{mviewFlattedCh}
    
    PREDICATOR(mviewFlattedCh)

    emit: 
    pepSeg = PREDICATOR.out.pepSeg
    msaEpitopeReport = PREDICATOR.out.msaEpitopeReport
}