#!/usr/bin/env python

import sys
import csv
import collections

from model.Segments import Segment

import click
from Bio import SeqIO


@click.command(help="Aim of this program is generate fasta file with peptide sequences")
@click.argument("seg", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument("pep", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--prefix", "-pf")
@click.option("--outdir", "-o", default="outdir")
def main(seg, pep, outdir:str, prefix:str):

    _setup_csv_field_size_limit()
    _read_mapped_segment_tsv(seg, pep, outdir, prefix)


def _setup_csv_field_size_limit():
    # increasing csv file size limit
    csv.field_size_limit(int(sys.maxsize/1000))

def _create_output_file_name(outdir, prefix, segment):
    return outdir + "/" + prefix + "_" +segment.get_name() + ".fasta"


def _read_mapped_segment_tsv(map_seg_file, pep, outdir, prefix, with_freq=False):
    DELIM="\t"
    SCF_POS=0 
    START_POS=1
    END_POS=2
    TOTAL_READS_MAPPED_POS=3
    MIN_MAPQ=4
    MAX_MAPQ=5
    AVG_MAPQ=6
    MEDIAN_MAPQ=7
    LIST_OF_READS=8

    segments = []

    with open(map_seg_file) as file:
        file_reader = csv.reader(file, delimiter=DELIM)
        for line in file_reader:
            s = Segment(line[SCF_POS], line[START_POS], line[END_POS], line[TOTAL_READS_MAPPED_POS] ,line[MIN_MAPQ], line[MAX_MAPQ], line[AVG_MAPQ], line[MEDIAN_MAPQ], line[LIST_OF_READS].split(";") )
            
            if with_freq:
                sequences_with_frequency = collections.Counter(s.get_peptide_ids())
                sequences = _build_list_of_sequences_with_frequency(dict(sequences_with_frequency), pep, True)
            else:
                distinct_sequences =  s.get_distinct_peptide_ids()
                sequences = _build_list_of_sequences(distinct_sequences, pep, True)

            output_name = _create_output_file_name(outdir, prefix, s)
            print(f"{len(sequences)} sequences was found")
            _write_output_file_(sequences, output_name)
            segments.append(s)
            # print(line[SCF_POS],line[START_POS], line[END_POS], line[TOTAL_READS_MAPPED_POS], line[AVG_MAPQ], line[MEDIAN_MAPQ], line[MIN_MAPQ], line[MAX_MAPQ], line[LIST_OF_READS])
            # bedtools merge  -c (columns) 1,5,5,5,5,1 -delim ";" -o count,min,max,mean,median,collapse > ${outFileName} 
    return segments

def _build_list_of_sequences_with_frequency(sequences_with_frequency, target_fasta_input, del_desc):
    record_dict = _read_fasta(target_fasta_input)
    sequences = []
    for seq_id, frequency in sequences_with_frequency.items():
        try:
            r = record_dict[seq_id]
            if del_desc:
                r.description = ''
            for _ in range(int(frequency)):
                sequences.append(r)
        except Exception as e:
            print("Ops", e, "occurred")
        
    return sequences

def _build_list_of_sequences(seq_id_list, target_fasta_input, rd):
    record_dict = _read_fasta(target_fasta_input)
    sequences = []
    for seq_id in seq_id_list:
        try:
            r = record_dict[seq_id]
            if rd:
                r.description = ''
            sequences.append(r)
        except Exception as e:
            print("Ops", e, "occurred")

    return sequences

def _read_fasta(faa_file):
    return SeqIO.index(faa_file, "fasta")

def _write_output_file_(sequences, output_file_name):
    with open(output_file_name, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

if __name__ == "__main__":
    main()
