include { PRE_EPITOPE_RETRIEVER } from './modules/local/pre_epitope_retriever/main.nf'

workflow preEpitopeRetriever{
    take:
    preEpitope
    
    main:
    PRE_EPITOPE_RETRIEVER(preEpitope.collect())
}
