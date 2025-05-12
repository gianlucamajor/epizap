
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
    
    def get_feature(self, scaffold, start, end, features): 
        annotation_found = []
        for annotation in self.segmentAnnotations:
            if (annotation.feature in features and annotation.get_segments_in_range(scaffold, start, end)):
                    target_segment_length = int(end - start)
                    segmentAnnotationLenth = int(annotation.end - annotation.start)
                    coverageAnnotation = float(segmentAnnotationLenth / target_segment_length)
                    if annotation.feature == "CDS":
                        product =  self._safe_unquote(annotation.get_attribute_by_name('product')) 
                        annotation_found.append(F"{annotation.locus} | {annotation.get_attribute_by_name('protein_id')} | {product} | {coverageAnnotation}")
                    if annotation.feature == "pseudogene":
                        note = self._safe_unquote(annotation.get_attribute_by_name('Note'))
                        annotation_found.append(F"{annotation.locus} | {annotation.get_attribute_by_name('gene_biotype')} | {note} | {coverageAnnotation}")
        return annotation_found
    
    def _safe_unquote(self, str_to_unquote):
        if str_to_unquote is not None:
            try:
                return unquote(str_to_unquote)
            except Exception as e:
                print(f"Error unquoting string: {str_to_unquote}, Error: {e}")
        else:
            return ""
