#!/usr/bin/env python

import sys
import csv

from model.Segment import Segment

import click
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing.pool import Pool


@click.command(help="Aim of this program is generate fasta file with peptide sequences")
@click.argument("seg", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument("pep", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--prefix", "-pf")
@click.option("--outdir", "-o", default="outdir")
@click.option("--threads", "-t", type=int, default=1)
def main(seg, pep, outdir:str, prefix:str, threads:int):
    _setup_csv_field_size_limit()

    clear_read_description = True
    segment_list = _read_mapped_segment_tsv(seg)
    pep_dic = _get_dic_from_fasta(pep)
    
    list_of_parameter =[
        (segment.get_distinct_peptide_ids(), prefix, segment.get_name(), outdir, pep_dic, clear_read_description) for segment in segment_list
    ]

    with Pool(threads) as pool:
        pool.starmap(_build_list_of_sequences_and_write_output, list_of_parameter)

 

def _setup_csv_field_size_limit():
    # increasing csv file size limit
    csv.field_size_limit(int(sys.maxsize/1000))

def _create_output_file_name(outdir, prefix, segment_name, total_of_sequence):
    return outdir + "/" + prefix + "_" + segment_name + "_" + str(total_of_sequence) + ".fasta"


def _read_mapped_segment_tsv(map_seg_file):
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
            segments.append(s)
            # print(line[SCF_POS],line[START_POS], line[END_POS], line[TOTAL_READS_MAPPED_POS], line[AVG_MAPQ], line[MEDIAN_MAPQ], line[MIN_MAPQ], line[MAX_MAPQ], line[LIST_OF_READS])
            # bedtools merge  -c (columns) 1,5,5,5,5,1 -delim ";" -o count,min,max,mean,median,collapse > ${outFileName}
    return segments
    

def _build_list_of_sequences_and_write_output(seq_id_list, prefix, segment_name, outdir, pep_record_dict, rd):    
    sequences = []
    for seq_id in seq_id_list:
        try:
            r = pep_record_dict[seq_id] 
            if rd:
                r.description = ''
            sequences.append(r)
        except Exception as e:
            print("Ops", e, "occurred")
    total_of_sequence = len(sequences)
    output_name = _create_output_file_name(outdir, prefix, segment_name, total_of_sequence)
    print(f"{total_of_sequence} sequences was found")
    _write_output_file_(sequences, output_name)

    return sequences

def _get_dic_from_fasta(file):
    p_dic = _read_fasta(file)
    p_keys =  p_dic.keys()
    new_dic = {}
    for k in p_keys:
        new_dic[k] = p_dic[k]

    return new_dic

def _read_fasta(faa_file):
    return SeqIO.index(faa_file, "fasta")

def _write_output_file_(sequences, output_file_name):
    with open(output_file_name, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

if __name__ == "__main__":
    main()
