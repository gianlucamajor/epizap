#!/usr/bin/env python

import os

import pysam
import click


SUPPLEMENTARY_SOCORE_TAG = "XS"
ALIGNMENT_SCORE_TAG = "AS"

@click.command(help="Aim of this program is extract all maximum score mappings from BAM file")
@click.argument("input", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument("output", type=click.Path(exists=True, file_okay=False, writable=True, dir_okay=True))
@click.option("--sufix", "-sfx", type=str, default="msm")
def main(input:click.Path, output:click.Path, sufix:str):    
    pysam.index(input)
    samfile = pysam.AlignmentFile(input, "rb")
    output_file_path =  get_output_file_path(input, output, sufix)
    new_samfile = pysam.AlignmentFile(output_file_path, "wb", template=samfile)

    primary_reads = []
    other_max_score_reads = []
    other_max_score_mapping_dict = {}
    for read in samfile.fetch():
        if not read.is_secondary:
            primary_reads.append(read)
            new_samfile.write(read)
            other_max_score_mapping = get_max_score_mappings(read)
            if  other_max_score_mapping is not None:
                other_max_score_mapping_dict[other_max_score_mapping['query_name']] = other_max_score_mapping['max_score_value']        

    for read in samfile.fetch():
        if read.is_secondary:
            if read.query_name in other_max_score_mapping_dict.keys():
                aligment_max_score = int(other_max_score_mapping_dict[read.query_name])
                read_alignment_score = int(read.get_tag(ALIGNMENT_SCORE_TAG))
                if read_alignment_score == aligment_max_score:
                    other_max_score_reads.append(read)
                    new_samfile.write(read)
    
    new_samfile.close()         
    samfile.close()

def get_output_file_path(input, outdir, suffix):
    input_file_name =  os.path.basename(input)
    input_file_name_splited =  os.path.splitext(input_file_name)
    output_file_base_name = input_file_name_splited[0]
    output_ext = input_file_name_splited[1]
    
    output_file_name =  f"{output_file_base_name}_{suffix}{output_ext}" 
    
    return os.path.join(outdir, output_file_name)
    
    

def get_max_score_mappings(read):
    if read.has_tag(SUPPLEMENTARY_SOCORE_TAG): 
        primary_score_value  = int(read.get_tag(ALIGNMENT_SCORE_TAG))
        supplementary_score_value = int(read.get_tag(SUPPLEMENTARY_SOCORE_TAG))        
        if supplementary_score_value == primary_score_value:
            return {'query_name':read.query_name, 'max_score_value': primary_score_value}


if __name__ == "__main__":
    main()
