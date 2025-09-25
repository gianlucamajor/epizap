class BlastResults:
    def __init__(self, blast_file):
        self.hits = []
        with open(blast_file) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.strip().split('\t')
                # blastp output format: -outfmt "6 std sallseqid qseq sseq qlen slen qcovs nident"
                # Adjust these indices if your file format is different
                # std: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
                # Here: 0-0	394449	57.143	7	3	0	2	8	1	7	31	13.1	394449	THRCRPA	THTSTPA	9	15	78	4
                self.hits.append(fields)

    def filter_hits(self, min_length=0, min_identity=0.0, no_gaps=False, qseqid=None):
        filtered = []
        for fields in self.hits:
            try:
                identity = float(fields[2])
                length = int(fields[3])
                gaps = int(fields[5])
            except (IndexError, ValueError):
                continue

            if length < min_length:
                continue
            if identity < min_identity:
                continue
            if no_gaps and gaps > 0:
                continue
            if qseqid and fields[0] != qseqid:
                continue

            filtered.append(fields)
        return filtered

    #WARNING: this function returns a best hit per sseqid, not per qseqid!
    def best_hits_by_sseqid(self, min_length=0, min_identity=0.0, no_gaps=False, qseqid=None):
        best_hits = {}
        for fields in self.filter_hits(min_length, min_identity, no_gaps, qseqid):
            sseqid = fields[1]
            identity = float(fields[2])
            length = int(fields[3])
            if sseqid not in best_hits:
                best_hits[sseqid] = fields
            else:
                prev_identity = float(best_hits[sseqid][2])
                prev_length = int(best_hits[sseqid][3])
                if (identity > prev_identity) or (identity == prev_identity and length > prev_length):
                    best_hits[sseqid] = fields
        return best_hits.values()