
include { MSA } from './modules/local/msa/main.nf'

workflow msa {
    take:
    pepSegCh

    main:
    pepSegCh
        .flatten()
        .filter({ it.countFasta() > 1 && it.countFasta() < 1000})
        .map{it ->
            def metaAlign = [id: it.name.replaceFirst(".fasta", ""), records: it.countFasta(), algorithm: "align"]
            [metaAlign, it]
        }
        .set{msaPepSegCh}    
    
    
    pepSegCh
        .flatten()
        .filter({ it.countFasta() >= 1000})
        .map{it ->
            def metaSuper5 = [id: it.name.replaceFirst(".fasta", ""), records: it.countFasta(), algorithm: "super5"]
            [metaSuper5, it]
        }
        .set{msaPepSegLargeCh}    
    
    msaPepSegCh.mix(msaPepSegLargeCh).set{allMSAPepSeg}
    
    
    MSA(allMSAPepSeg) 

    emit: 
    MSA.out

}

