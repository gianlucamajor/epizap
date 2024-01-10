#!/usr/bin/env python

import sys
import csv
import collections

from model.Segment import Segment

import click
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@click.command(help="Aim of this program is generate fasta file with peptide sequences")
@click.argument("seg", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument("pep", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--prefix", "-pf")
@click.option("--outdir", "-o", default="outdir")
@click.option("--with_frequency", "-with_freq", default=False)
def main(seg, pep, outdir:str, prefix:str, with_frequency:bool):

    _setup_csv_field_size_limit()
    _read_mapped_segment_tsv(seg, pep, outdir, prefix, with_frequency)


def _setup_csv_field_size_limit():
    # increasing csv file size limit
    csv.field_size_limit(int(sys.maxsize/1000))

def _create_output_file_name(outdir, prefix, segment, total_of_sequence):
    return outdir + "/" + prefix + "_" +segment.get_name()+ "_" + str(total_of_sequence) + ".fasta"


def _read_mapped_segment_tsv(map_seg_file, pep, outdir, prefix, with_frequency):
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

    with open(map_seg_file) as file:
        file_reader = csv.reader(file, delimiter=DELIM)
        for line in file_reader:
            s = Segment(line[SCF_POS], line[START_POS], line[END_POS], line[TOTAL_READS_MAPPED_POS] ,line[MIN_MAPQ], line[MAX_MAPQ], line[AVG_MAPQ], line[MEDIAN_MAPQ], line[LIST_OF_READS].split(";") )
            
            pep_record_dict = _read_fasta(pep)
            
            if with_frequency:
                sequences_with_frequency = collections.Counter(s.get_peptide_ids())
                sequences = _build_list_of_sequences_with_frequency(dict(sequences_with_frequency), pep_record_dict, True)
            else:
                distinct_sequences =  s.get_distinct_peptide_ids()
                sequences = _build_list_of_sequences(distinct_sequences, pep_record_dict, True)

            total_of_sequence = len(sequences)
            output_name = _create_output_file_name(outdir, prefix, s, total_of_sequence)
            print(f"{total_of_sequence} sequences was found")
            _write_output_file_(sequences, output_name)
            
            # print(line[SCF_POS],line[START_POS], line[END_POS], line[TOTAL_READS_MAPPED_POS], line[AVG_MAPQ], line[MEDIAN_MAPQ], line[MIN_MAPQ], line[MAX_MAPQ], line[LIST_OF_READS])
            # bedtools merge  -c (columns) 1,5,5,5,5,1 -delim ";" -o count,min,max,mean,median,collapse > ${outFileName}
    

def _build_list_of_sequences_with_frequency(sequences_with_frequency, pep_record_dict, del_desc):
    sequences = []
    for seq_id, frequency in sequences_with_frequency.items():
        try:
            r = pep_record_dict[seq_id]
            if del_desc:
                r.description = ''
            for idx in range(int(frequency)):
                new_sequence_id = r.id
                if idx > 0:
                    new_sequence_id = seq_id + f"_clone_{idx}"
                new_record = SeqRecord(
                    r.seq,
                    id=new_sequence_id,
                    name=r.name,
                    description=r.description)

                sequences.append(new_record)
        except Exception as e:
            print("Ops", e, "occurred")
    return sequences

def _build_list_of_sequences(seq_id_list, pep_record_dict, rd):
    sequences = []
    for seq_id in seq_id_list:
        try:
            r = pep_record_dict[seq_id]
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
