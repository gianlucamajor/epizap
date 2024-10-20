import pysam

bFilePath = "results/mapping/CCC_severe_a.bam"

pysam.index(bFilePath)
samfile = pysam.AlignmentFile(bFilePath, "rb")

for read in samfile.fetch():
    if not read.is_secondary:
        print(read.get_tags())
        # print(read)

samfile.close()