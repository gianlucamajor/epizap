#!/usr/bin/env python

import os
import re
import csv

import click
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# python3 bin/epitope_predictor.py results_05_11_2024/mview/*.html -o results_05_11_2024/

CONSENSUS_PATTERN = re.compile(r'^consensus/(\d{2,3})%')
EPITOPE_PATTERN = re.compile(r'[A-Z]{8,}')
# IDENTITY_AND_COVERAGE_PATTERN = re.compile(r'^\d+\s+\d+\s+\d{1,3}\.\d%?\s+\d{1,3}\.\d%?') # 
IDENTITY_AND_COVERAGE_PATTERN = re.compile(r'^\d+\s+\d+-?\d*\s+\d{1,3}\.\d%?\s+\d{1,3}\.\d%?') # with hyphen
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
        antigen_consensus = get_antigen_consensus(msa_epitopes['epitope_candidates'])
        # verify_epitope(msa_epitopes, ee, antigen_consensus)

        msa_epitopes['antigen_conseus'] = antigen_consensus
        msa_epitopes['epitope_elected'] = ee
        msa_epitopes_list.append(msa_epitopes)

    msa_epitopes_list_sorted = sorted(msa_epitopes_list, key=lambda epitope: epitope['msa_avg_identity'], reverse=True)
    epitope_writer(msa_epitopes_list_sorted, outdir)


def epitope_writer(epitope_list, outdir):
    tsv_lines_report = []
    epitopes_seq_list = []
    antigen_consesus_list = []
 
    for ept in epitope_list:
        name = ept['file_name'].split(".html")[0] # to remove .html
        number_of_epitopes = 0
        consensus = None
        epitopes = None
        epitopes_by_cc = []
        
        if ept['epitope_elected']:
            number_of_epitopes = len(ept['epitope_elected']['epitopes'])
            consensus = ept['epitope_elected']['consensus']
            epitopes = ept['epitope_elected']['epitopes']
            antigen_consesus = ept['antigen_conseus']
            antigen_consesus_list.append(_create_seq_record_with_name(antigen_consesus, name))
            for idx, e in enumerate(epitopes):
                seq_rec = _create_seq_record(e, name, idx)
                epitopes_by_cc.append(_create_simple_sequence(e, name, idx))
                epitopes_seq_list.append(seq_rec)
                
                
        # "cc_id", "number_of_epitopes", "epitopes", "consensus_percentage", "msa_avg_identity"
        tsv_lines_report.append(f"{_remove_cc_acronym(name)}\t{number_of_epitopes}\t{epitopes_by_cc}\t{consensus}\t{ept['msa_avg_identity']}")

        tsv_of = f"{outdir}/epitopes-msa-predict-report.tsv"
        fa_of = f"{outdir}/epitopes.fasta"
        antigen_of = f"{outdir}/antigens.fasta"

        _write_tsv_file(tsv_of, tsv_lines_report)
        _write_output_fasta_file(epitopes_seq_list, fa_of)
        _write_output_fasta_file(antigen_consesus_list, antigen_of)

        


def _create_seq_record(sequence, name, idx):
    name_with_idx = _create_seq_name_with_idx(name, idx)
    return SeqRecord(Seq(sequence), id=name_with_idx, name=name_with_idx, description="")

def _create_seq_record_with_name(sequence, name):
    return SeqRecord(Seq(sequence), id=name, name=name, description="")

def _create_seq_name_with_idx(name, idx):
    name = _remove_cc_acronym(name) # to remove cc acronym
    name_with_idx = f"{name}-{idx}"
    return name_with_idx


def _create_simple_sequence(sequence, name, idx):
    name_with_idx = _create_seq_name_with_idx(name, idx)
    return {
    "id": name_with_idx,
    "name": name_with_idx,
    "seq": sequence }
    


def _write_tsv_file(tsv_of, tsv_lines_list):
    with open(tsv_of, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(["cc_id", "number_of_epitopes", "epitopes", "consensus_percentage", "msa_avg_identity"])
        for line in tsv_lines_list:
            tsv_line = line.split("\t")
            writer.writerow(tsv_line)
    
    

def _write_output_fasta_file(sequences, output_file_name):
    with open(output_file_name, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

def epitope_voter(msa_epitopes_candidate):
    if len(msa_epitopes_candidate) < 1:
        return None
    msa_epitopes_candidate_sorted = sorted(msa_epitopes_candidate, key=lambda epitope: epitope['consensus'], reverse=True)
    return msa_epitopes_candidate_sorted[0]

def get_antigen_consensus(msa_epitopes_candidate):
    if len(msa_epitopes_candidate) < 1:
        return None
    
    msa_epitopes_candidate_sorted = sorted(msa_epitopes_candidate, key=lambda epitope: epitope['consensus'], reverse=True)
    
    antigen_consensus = None
    for epitopes_and_consensus in msa_epitopes_candidate_sorted:
         epitopes_and_consensus_sorted_by_longest = sorted(epitopes_and_consensus['epitopes'], key=lambda epitope: len(epitope), reverse=True)
         longest_sequence = epitopes_and_consensus_sorted_by_longest[0]
         if antigen_consensus is None or len(longest_sequence) > len(antigen_consensus):
                antigen_consensus = longest_sequence
                
    return antigen_consensus
def verify_epitope(msa_epitopes, ee, antigen_consensus):
    

    if ee is not None:
        for epitope_sequence in ee['epitopes']:
            if epitope_sequence not in antigen_consensus:
                print(f"{msa_epitopes['file_name']} {msa_epitopes['msa_avg_identity']}")
                print(f"{len(ee['epitopes'])} - {ee['epitopes']}")
                print(f"WARNING: {epitope_sequence}  is not in the antigen consensus {antigen_consensus}")


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

def _remove_cc_acronym(name):
    # to remove cc acronym
    if name.startswith("cc-"):
        name = name.split("-")[1]
    return name

if __name__ == "__main__":
    main()