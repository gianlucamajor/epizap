#!/usr/bin/env python

class SegmentAnnotation:
    def __init__(self, scaffold, start, end, reads_count, min_mapq, reads_strand, max_mapq, avg_mapq, median_mapq, scaffold2, source, feature, start2, end2, score, annotation_strand, phase, attributes):
        self.locus = f"{scaffold}-{start}-{end}"
        self.scaffold = scaffold
        self.start = int(start)
        self.end = int(end)
        self.reads_count = int(reads_count)
        self.min_mapq = int(min_mapq)
        self.reads_strand = reads_strand
        self.max_mapq = int(max_mapq)
        self.avg_mapq = float(avg_mapq)
        self.median_mapq = float(median_mapq)
        self.scaffold2 = scaffold2
        self.source = source
        self.feature = feature
        self.start2 = int(start2)
        self.end2 = int(end2)
        self.score = score
        self.annotation_strand = annotation_strand
        self.phase = phase
        self.attributes = attributes
    
    def get_attribute_by_name(self, name):
        attributes_dict = dict(attribute.split('=') for attribute in self.attributes.split(';') if '=' in attribute)
        return attributes_dict.get(name)

    def get_attributes(self):
        return self.attributes
    
    def get_segments_in_range(self, scaffold, start, end):
        if self.scaffold == scaffold and self.start >= start and self.end <= end:
            return self
        return None

    def __str__(self):
        return f"{self.scaffold}\t{self.start}\t{self.end}\t{self.reads_count}\t{self.min_mapq}\t{self.max_mapq}\t{self.avg_mapq}\t{self.median_mapq}\t{self.reads_strand}\t{self.scaffold2}\t{self.source}\t{self.feature}\t{self.start2}\t{self.end2}\t{self.score}\t{self.annotation_strand}\t{self.phase}\t{self.attributes}"
