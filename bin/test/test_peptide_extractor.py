import unittest
import os
import sys
import csv 

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.append("bin")

import peptide_extractor as main


class TestPeptideExtractor(unittest.TestCase):

    def test_setup_csv_field_size_limit(self):
        main._setup_csv_field_size_limit()
        self.assertEqual(csv.field_size_limit(), int(sys.maxsize/1000))
    
    def test_build_list_of_sequences(self):
        # Arrange
        seq_id_list = ["seq1", "seq2", "seq3"]
        pep_record_dict = {"seq1": SeqRecord(Seq("ACGT"), id="seq1", name="seq1", description="This is a test sequence."),
                          "seq2": SeqRecord(Seq("ACGT"), id="seq2", name="seq2", description="This is another test sequence."),
                          "seq3": SeqRecord(Seq("ACGT"), id="seq3", name="seq3", description="This is a third test sequence.")}
        rd = True

        # Act
        sequences = main._build_list_of_sequences(seq_id_list, pep_record_dict, rd)

        # Assert
        self.assertEqual(len(sequences), 3)
        self.assertEqual(sequences[0].id, "seq1")
        self.assertEqual(sequences[0].seq, Seq("ACGT"))
        self.assertEqual(sequences[0].description, "")
        self.assertEqual(sequences[1].id, "seq2")
        self.assertEqual(sequences[1].seq, Seq("ACGT"))
        self.assertEqual(sequences[1].description, "")
        self.assertEqual(sequences[2].id, "seq3")
        self.assertEqual(sequences[2].seq, Seq("ACGT"))
        self.assertEqual(sequences[2].description, "")


if __name__ == "__main__":
    unittest.main()