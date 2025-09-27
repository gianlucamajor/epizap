from Bio import AlignIO
# from Bio.Align import MultipleSeqAlignment
# from Bio.Align.AlignInfo import SummaryInfo
# Assuming your alignment file is named "aligned_sequences.fasta"
from Bio.Align import AlignInfo

msa = AlignIO.read("/home/gianluca/workspace/epizap/results_30_04_2025/msa/cc-1391-msa.msc", "fasta")
summary_align = AlignInfo.SummaryInfo(msa)

alignment = msa.alignment
# Generate consensus with a threshold of 0.7 (70% frequency)
consensus_sequence = summary_align.dumb_consensus(0.501, "X")


print(alignment)
print(consensus_sequence)
