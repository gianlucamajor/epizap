#!/usr/bin/env python

import os
import re

import click
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# python3 bin/epitope_predictor.py results_05_11_2024/mview/*.html -o results_05_11_2024/

CONSENSUS_PATTERN = re.compile(r'^consensus/(\d{2,3})%')
EPITOPE_PATTERN = re.compile(r'[A-Z]{5,}')
IDENTITY_AND_COVERAGE_PATTERN = re.compile(r'^\d+\s+\d+\s+\d{1,3}\.\d%?\s+\d{1,3}\.\d%?')
PERCENTAGE_OF_IDENTITY_PATTERN = re.compile(r'(\d{1,3}\.\d)(?=%$)')

@click.command(help="Aim of this program is predict epitope from the multiple sequence align (MSA) core.")
@click.argument("msa_files", nargs=-1, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--outdir", "-o", type=click.Path(), default="")
def main(msa_files:click.Path, outdir:click.Path):
    

    msa_epitopes_list = []
    for msa_file in msa_files:
        msa_epitopes = to_process_msa(msa_file)
        
        # epitope voter 
        ee =  epitope_voter(msa_epitopes['epitope_candidates'])
        msa_epitopes['epitope_elected'] = ee
        msa_epitopes_list.append(msa_epitopes)

    msa_epitopes_list_sorted = sorted(msa_epitopes_list, key=lambda epitope: epitope['msa_avg_identity'], reverse=True)
    epitope_writer(msa_epitopes_list_sorted, outdir)




def epitope_writer(epitope_list, outdir):
    tsv_lines_report = []
    epitopes_seq_list = []
 
    for ept in epitope_list:
        name = ept['file_name'].split(".html")[0] # to remove .html
        number_of_epitopes = 0
        consensus = None
        epitopes = None
        
        if ept['epitope_elected']:
            number_of_epitopes = len(ept['epitope_elected']['epitopes'])
            consensus = ept['epitope_elected']['consensus']
            epitopes = ept['epitope_elected']['epitopes']
            
            for e in epitopes:
                seq_rec = _create_seq_record(e, name, epitopes.index(e))
                epitopes_seq_list.append(seq_rec)
                
                

        tsv_lines_report.append(f" {name} \t {number_of_epitopes} \t {epitopes} \t {consensus} \t {ept['msa_avg_identity']}")
    
        of = f"{outdir}/epitopes.fasta"
        _write_output_file_(epitopes_seq_list, of)

        


def _create_seq_record(sequence, name, idx):
    name_with_idx = f"{name}-{idx}"
    return SeqRecord(Seq(sequence), id=name_with_idx, name=name_with_idx, description="")


def _write_output_file_(sequences, output_file_name):
    with open(output_file_name, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

def epitope_voter(msa_epitopes_candidate):
    if len(msa_epitopes_candidate) < 1:
        return None
    msa_epitopes_candidate_sorted = sorted(msa_epitopes_candidate, key=lambda epitope: epitope['consensus'], reverse=True)
    return msa_epitopes_candidate_sorted[0]


def to_process_msa(msa_file):
    with open(msa_file, 'r') as msa_file_opened:
        msa_epitopes = {}
        epitopes_candidate_list = []
        identity_and_coverage_list = []

        for msa_line in msa_file_opened:
            msa_line_striped = msa_line.strip()

            identity_and_coverage_found = IDENTITY_AND_COVERAGE_PATTERN.search(msa_line_striped)
            if identity_and_coverage_found:
                identity_and_coverage_list.append(identity_and_coverage_found.group(0))

            consensus_found = CONSENSUS_PATTERN.search(msa_line_striped)
            if consensus_found:
                consensus_value = consensus_found.group(1)                
                epitopes_found = EPITOPE_PATTERN.findall(msa_line_striped)
                if epitopes_found:
                    epitopes_candidate = {}
                    epitopes_candidate['consensus'] = int(consensus_value)
                    epitopes_candidate['epitopes'] = epitopes_found

                    epitopes_candidate_list.append(epitopes_candidate)

        msa_epitopes['file_name'] = os.path.basename(msa_file_opened.name)
        msa_epitopes['epitope_candidates'] = epitopes_candidate_list
        msa_epitopes['msa_avg_identity'] = average_identity_calculator(identity_and_coverage_list)
        return msa_epitopes

def average_identity_calculator(identity_and_coverage_list):
    size = len(identity_and_coverage_list)
    total = 0
    for sequence_line in identity_and_coverage_list:
        pid_found = PERCENTAGE_OF_IDENTITY_PATTERN.search(sequence_line)
        pid_value = pid_found.group(0)
        total += float(pid_value)
    avg_identity = total/size

    return avg_identity
    
                

                


if __name__ == "__main__":
    main()