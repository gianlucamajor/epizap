import pysam

bFilePath = "results/mapping/asympto_a.bam"
# bFilePath = "results_19_10_2024/mapping/control_and_chagasic_patients.bam"


pysam.index(bFilePath)
samfile = pysam.AlignmentFile(bFilePath, "rb")
new_samfile = pysam.AlignmentFile("new_asympto_a.bam", "wb", template=samfile)
primary_reads = []
for read in samfile.fetch():
    if not read.is_secondary:
        primary_reads.append(read)
        new_samfile.write(read)

new_samfile.close()    
samfile.close()

print(len(primary_reads))
for pr in primary_reads:

    print(pr.get_tag("AS"), pr.get_tag("XS"), pr.get_tags())

