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
    PREDICATOR.out
}