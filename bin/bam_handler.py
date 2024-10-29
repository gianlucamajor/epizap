import pysam



# bFilePath = "results_19_10_2024/mapping/control_and_chagasic_patients.bam"
# bNewFilePath = "results_19_10_2024/mapping/new_control_and_chagasic_patients.bam"

bFilePath = "results/mapping/asympto_a.bam"
bNewFilePath = "results/mapping/new_asympto_a.bam"


pysam.index(bFilePath)
samfile = pysam.AlignmentFile(bFilePath, "rb")
new_samfile = pysam.AlignmentFile(bNewFilePath, "wb", template=samfile)


SUPPLEMENTARY_SOCORE_TAG = "XS"
ALIGNMENT_SCORE_TAG = "AS"

def get_max_score_mappings(read):
    if read.has_tag(SUPPLEMENTARY_SOCORE_TAG): 
        primary_score_value  = int(read.get_tag(ALIGNMENT_SCORE_TAG))
        supplementary_score_value = int(read.get_tag(SUPPLEMENTARY_SOCORE_TAG))        
        if supplementary_score_value == primary_score_value:
            return {'query_name':read.query_name, 'max_score_value': primary_score_value}

primary_reads = []
other_max_score_reads = []
other_max_score_mapping_dict = {}
for read in samfile.fetch():
    if not read.is_secondary:
        primary_reads.append(read)
        new_samfile.write(read)
        other_max_score_mapping = get_max_score_mappings(read)
        if  other_max_score_mapping is not None:
            other_max_score_mapping_dict[other_max_score_mapping['query_name']] = other_max_score_mapping['max_score_value']        

for read in samfile.fetch():
    if read.is_secondary:
        if read.query_name in other_max_score_mapping_dict.keys():
            aligment_max_score = int(other_max_score_mapping_dict[read.query_name])
            read_alignment_score = int(read.get_tag(ALIGNMENT_SCORE_TAG))
            if read_alignment_score == aligment_max_score:
                other_max_score_reads.append(read)
                new_samfile.write(read)


new_samfile.close()         
samfile.close()



print(len(primary_reads))
print(len(other_max_score_reads))
print(len(other_max_score_reads) + len(primary_reads))

