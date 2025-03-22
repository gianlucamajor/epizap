
from model.SegmentAnnotation import SegmentAnnotation

from urllib.parse import unquote

class AnnotationHandler:
    def __init__(self, file_path):
        self.file_path = file_path
        self.segmentAnnotations = self._load_file()

    def _load_file(self):
         segmentAnnotations = []
         with open(self.file_path, 'r') as file:
            for line in file:
                columns = line.strip().split('\t')
                segmentAnnotations.append(SegmentAnnotation(*columns))
         return segmentAnnotations
    
    def get_feature(self, scaffold, start, end, feature): 
        annotation_found = []
        for annotation in self.segmentAnnotations:
            if (annotation.feature == feature and annotation.get_segments_in_range(scaffold, start, end)):
                    target_segment_length = int(end - start)
                    segmentAnnotationLenth = int(annotation.end - annotation.start)
                    coverageAnnotation = float(segmentAnnotationLenth / target_segment_length)
                    annotation_found.append(F"{annotation.locus} | {annotation.get_attribute_by_name('protein_id')} | {unquote(annotation.get_attribute_by_name('product'))} | {coverageAnnotation}")
        return annotation_found

    def print_annotations(self):
        for annotation in self.segmentAnnotations:
            print(annotation)

# if __name__ == "__main__":
#     annHandler = AnnotationHandler("/home/gianluca/workspace/epizap/results_graph_21_01_2025/segments-annotated/control_and_chagasic_patients_mapped-segment-annotated.tsv")
#     # annHandler.print_annotations()
#     annHandler.get_feature("CM026583.1", 913586, 914466, "CDS")