
include { MSA } from './modules/local/msa/main.nf'
include { HMM_PROFILE_BUILDER } from './modules/local/hmm_profile_builder/main.nf'
include { HMM_CONSENSUS_BUILDER } from './modules/local/hmm_consensus_builder/main.nf'

workflow consensusBuilder {
    take:
    pepSegCh

    main:
    pepSegCh
        .flatten()
        .filter({ it.countFasta() > 1 && it.countFasta() < 1000})
        .map{it ->
            metaAlign = [id: it.name.replaceFirst(".fasta", ""), records: it.countFasta(), algorithm: "align"]
            [metaAlign, it]
        }
        .set{msaPepSegCh}    
    
    
    pepSegCh
        .flatten()
        .filter({ it.countFasta() >= 1000})
        .map{it ->
            metaSuper5 = [id: it.name.replaceFirst(".fasta", ""), records: it.countFasta(), algorithm: "super5"]
            [metaSuper5, it]
        }
        .set{msaPepSegLargeCh}    
    
    msaPepSegCh.mix(msaPepSegLargeCh).set{allMSAPepSeg}
    
    
    MSA(allMSAPepSeg) | HMM_PROFILE_BUILDER | HMM_CONSENSUS_BUILDER

    emit: 
    HMM_CONSENSUS_BUILDER.out

}

