include { JOIN_MSA_CORE_AND_LONELY_EPITOPES } from './modules/local/join_msa_core_and_lonely_epitopes/main.nf'

workflow msaCoreAndLonelyEpitopesJoiner{
    take:
    msaCoreAndLonelyEpitopes
    
    main:
    JOIN_MSA_CORE_AND_LONELY_EPITOPES(msaCoreAndLonelyEpitopes.collect())
}
