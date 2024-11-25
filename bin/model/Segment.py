#!/usr/bin/env python

class Segment():
    
    def __init__(self, scaffod:str, start:int, end:int, count:int, min_mapq:float, max_mapq:float, avg_mapq:float, median_mapq:float, reads:[]) -> None:
        self.scaffod = scaffod
        self.start = start
        self.end = end
        self.reads_count = count
        self.min_mapq = min_mapq
        self.max_mapq = max_mapq
        self.avg_mapq = avg_mapq
        self.median_mapq = median_mapq
        self.reads = reads

    def __str__(self) -> str:
        return f"{self.scaffod}-{self.start}-{self.end}-{len(self.reads)}"
    
    def get_name(self) -> str:
        return f"{self.scaffod}-{self.start}-{self.end}-{len(self.reads)}"
    
    def get_peptide_ids(self) -> list:
        return [r.split('_')[1] for r in self.reads]
    
    def get_distinct_peptide_ids(self) -> set:
        return set(self.get_peptide_ids())
    
    def get_set_of_reads(self) -> set:
        return set(self.reads)
    


