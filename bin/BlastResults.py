from enum import IntEnum

class BlastColumns(IntEnum):
    """Enum for BLAST output column indices"""
    QSEQID = 0    # Query sequence id
    SSEQID = 1    # Subject sequence id
    PIDENT = 2    # Percentage identity
    LENGTH = 3    # Alignment length
    MISMATCH = 4  # Number of mismatches
    GAPOPEN = 5   # Number of gap openings
    QSTART = 6    # Start of alignment in query
    QEND = 7      # End of alignment in query
    SSTART = 8    # Start of alignment in subject
    SEND = 9      # End of alignment in subject
    EVALUE = 10   # Expect value
    BITSCORE = 11 # Bit score
    QSEQ = 13     # Query sequence
    SSEQ = 14     # Subject sequence
    STITLE = 19   # Subject title 

class BlastResults:
    def __init__(self, blast_file):
        self.hits = []
        with open(blast_file) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.strip().split('\t')
                # blastp output format: -outfmt "6 std sallseqid qseq sseq qlen slen qcovs nident stitle"
                # Adjust these indices if your file format is different
                # std: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
                # Here: 10-0 PWV10101.1   100.000 10   0    0    1    10   152  161  8.7  22.3 PWV10101.1   QGAGNLSSGA   QGAGNLSSGA   10   166  100  10   PWV10101.1 putative RNase H [Trypanosoma cruzi]
                # The stile is only column: PWV10101.1 putative RNase H [Trypanosoma cruzi].
                self.hits.append(fields)

    def filter_hits(self, min_length=0, min_identity=0.0, no_gaps=False, qseqid=None):
        filtered = []
        for fields in self.hits:
            try:
                identity = float(fields[BlastColumns.PIDENT])
                length = int(fields[BlastColumns.LENGTH])
                gaps = int(fields[BlastColumns.GAPOPEN])
            except (IndexError, ValueError):
                continue

            if length < min_length:
                continue
            if identity < min_identity:
                continue
            if no_gaps and gaps > 0:
                continue
            if qseqid and fields[BlastColumns.QSEQID] != qseqid:
                continue

            filtered.append(fields)
        return filtered

    #WARNING: this function returns a best hit per sseqid, not per qseqid!
    def best_hits_by_sseqid(self, min_length=0, min_identity=0.0, no_gaps=False, qseqid=None):
        best_hits = {}
        for fields in self.filter_hits(min_length, min_identity, no_gaps, qseqid):
            sseqid = fields[BlastColumns.SSEQID]
            identity = float(fields[BlastColumns.PIDENT])
            length = int(fields[BlastColumns.LENGTH])
            if sseqid not in best_hits:
                best_hits[sseqid] = fields
            else:
                prev_identity = float(best_hits[sseqid][BlastColumns.PIDENT])
                prev_length = int(best_hits[sseqid][BlastColumns.LENGTH])
                if (identity > prev_identity) or (identity == prev_identity and length > prev_length):
                    best_hits[sseqid] = fields
        return best_hits.values()

    def best_hits_by_qseqid_and_bitscore(self, min_length=0, min_identity=0.0, no_gaps=False, qseqid=None):
        best_hits = {}
        for fields in self.filter_hits(min_length, min_identity, no_gaps, qseqid):
            qseqid = fields[BlastColumns.QSEQID]
            bitscore = float(fields[BlastColumns.BITSCORE])
            if qseqid not in best_hits:
                best_hits[qseqid] = fields
            else:
                prev_bitscore = float(best_hits[qseqid][BlastColumns.BITSCORE])
                if bitscore > prev_bitscore:
                    best_hits[qseqid] = fields
                
        return best_hits
    