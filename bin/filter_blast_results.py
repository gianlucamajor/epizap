#!/usr/bin/env python

import click
from BlastResults import BlastResults
from IEDBEpitopeTableHandler import IEDBEpitopeTableHandler

# Define column names for clarity
QUERY_ID = 0
SUBJECT_ID = 1
IDENTITY = 2
LENGTH = 3
GAPS = 5
QUERY_SEQ = 13
SUBJECT_SEQ = 14

# THIS IS MOKEY CODE TO TEST IEDB TABLE HANDLER

@click.command(help="Filter BLAST tabular results by hit length, identity, gap presence, and optionally by qseqid. Returns only the best hit per sseqid.")
@click.argument("blast_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--min-length", type=int, default=0, help="Minimum length of the hit (alignment length).")
@click.option("--min-identity", type=float, default=0.0, help="Minimum percent identity.")
@click.option("--no-gaps", is_flag=True, help="Exclude hits with gaps in the alignment.")
@click.option("--qseqid", type=str, default=None, help="Filter by specific qseqid value.")
@click.option("--print-header", is_flag=True, help="Print header line.")
def main(blast_file, min_length, min_identity, no_gaps, qseqid, print_header):
    blast = BlastResults(blast_file)
    hits = blast.filter_hits(min_length, min_identity, no_gaps, qseqid)
    # iedbTableHandler = IEDBEpitopeTableHandler("dataSet/iedb/epitope_table_export_1747317719_15052025.csv")
    iedbTableHandler = IEDBEpitopeTableHandler("dataSet/iedb/human_epitope_table_export_1757095864.csv")

    
    if print_header:
        print("QueryID\tSubjectID\tIdentity\tLength\tGaps\tQuerySeq\tSubjectSeq")

    for fields in hits:
        #if result is not None:
        #     print(result['Epitope_id'])
        #     print(result['epitope_sequence'])
        #     print(result['Epitope - Source Molecule'])
        #     print(result['Epitope - Source Molecule IRI'])
        iedbEpitopeInfo = iedbTableHandler.get_by_epitope_id(fields[SUBJECT_ID])
        
        print(f"{fields[QUERY_ID]}\t{fields[SUBJECT_ID]}\t{fields[IDENTITY]}\t{fields[LENGTH]}\t{fields[GAPS]}\t{fields[QUERY_SEQ]}\t{fields[SUBJECT_SEQ]}\t{iedbEpitopeInfo['Epitope_id']}\t{iedbEpitopeInfo['epitope_sequence']}\t{iedbEpitopeInfo['Epitope - Source Molecule']}\t{iedbEpitopeInfo['Epitope - Source Molecule IRI']}")
        
if __name__ == "__main__":
    main()