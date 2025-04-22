class BlastResult:
    # -outfmt '6 std sallseqid qseq sseq qlen slen qcovs nident'
    def __init__(self, line: str):
        # Split the line into fields using tab as the delimiter
        fields = line.strip().split("\t")
        
        # Assign fields to attributes based on the provided column names
        self.qseqid = fields[0]       # Query sequence ID
        self.sseqid = fields[1]       # Subject sequence ID
        self.pident = float(fields[2])  # Percentage of identical matches
        self.length = int(fields[3])  # Alignment length
        self.mismatch = int(fields[4])  # Number of mismatches
        self.gapopen = int(fields[5])  # Number of gap openings
        self.qstart = int(fields[6])  # Start of alignment in query
        self.qend = int(fields[7])    # End of alignment in query
        self.sstart = int(fields[8])  # Start of alignment in subject
        self.send = int(fields[9])    # End of alignment in subject
        self.evalue = float(fields[10])  # Expect value
        self.bitscore = float(fields[11])  # Bit score
        self.sallseqid = fields[12]   # All subject sequence IDs
        self.qseq = fields[13]        # Aligned part of query sequence
        self.sseq = fields[14]        # Aligned part of subject sequence
        self.qlen = int(fields[15])   # Query sequence length
        self.slen = int(fields[16])   # Subject sequence length
        self.qcovs = int(fields[17])  # Query coverage per subject
        self.nident = int(fields[18]) # Number of identical matches

 

    def __str__(self):
        return (
            f"Line(qseqid={self.qseqid}, sseqid={self.sseqid}, pident={self.pident}, "
            f"length={self.length}, mismatch={self.mismatch}, gapopen={self.gapopen}, "
            f"qstart={self.qstart}, qend={self.qend}, sstart={self.sstart}, send={self.send}, "
            f"evalue={self.evalue}, bitscore={self.bitscore}, sallseqid={self.sallseqid}, "
            f"qseq={self.qseq}, sseq={self.sseq}, qlen={self.qlen}, slen={self.slen}, "
            f"qcovs={self.qcovs}, nident={self.nident})"
        )