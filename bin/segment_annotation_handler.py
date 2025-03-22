#!/usr/bin/env python

import click
from urllib.parse import unquote
from model.SegmentAnnotation import SegmentAnnotation

@click.command()
@click.argument('file_path', type=click.Path(exists=True))
@click.option('--scaffold', required=True, help='Scaffold to filter')
@click.option('--start', required=True, type=int, help='Start position to filter')
@click.option('--end', required=True, type=int, help='End position to filter')
@click.option('--feature', required=True, help='Feature to filter')
def main(file_path, scaffold, start, end, feature):
    read_tsv(file_path, scaffold, start, end, feature)

def read_tsv(file_path, scaffold, start, end, feature):
    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            columns = line.strip().split('\t')
            annotation = SegmentAnnotation(*columns)
            if (annotation.feature == feature and annotation.get_segments_in_range(scaffold, start, end)):
                target_segment_length = int(end - start)
                segmentAnnotationLenth = int(annotation.end - annotation.start)
                coverageAnnotation = float(segmentAnnotationLenth / target_segment_length)
                print(annotation.locus, annotation.get_attribute_by_name('protein_id'), unquote(annotation.get_attribute_by_name('product')), coverageAnnotation)

if __name__ == "__main__":
    main()
