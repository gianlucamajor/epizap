include { LONELY } from './modules/local/lonely/main.nf'

workflow lonelyPeptideRetriever{
    take:
    lonelyPepSegCh
    
    main:
    lonelyPepSegCh
        .flatten()
        .filter({ it.countFasta() < 2})
        .collect()
        .set{lonelyPepSegCh}
    
    LONELY(lonelyPepSegCh)

    emit: 
    LONELY.out
}
