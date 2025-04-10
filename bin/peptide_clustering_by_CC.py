#!/usr/bin/env python

import os

import click
import networkx as nx
import matplotlib.pyplot as plt
import pickle

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import logging

logger = logging.getLogger("PCG")
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


@click.command(help="Aim of this program is to identify the connected components of a graph and save the peptides of each connected component in a fasta file.")
@click.argument("graph_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument("peptide_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("output_dir", "-o", type=click.Path(exists=True, file_okay=False, dir_okay=True), default=None, help="The output directory to save the graph.")
@click.option("cc_output", "-cc-output", type=click.Path(exists=False, file_okay=False, dir_okay=True), default="connected-components", help="The output directory to save the peptides on connected components.")
def main(graph_file: click.Path, peptide_file: click.Path, output_dir: click.Path, cc_output: click.Path):
    
    output_path = _get_output_path(graph_file, output_dir)
    
    cc_out_dir = os.path.join(output_path, cc_output)
    os.makedirs(cc_out_dir, exist_ok=True)
    
    input_peptides = _read_fasta(peptide_file)


    with open(graph_file, 'rb') as f:
        graph = pickle.load(f)
    
    # Print basic information about the graph
    logger.info(f"Number of nodes: {graph.number_of_nodes()}")
    logger.info(f"Number of edges: {graph.number_of_edges()}")
    
    cc_with_only_one_peptide_list = []
    
    cc_list = identifying_and_sort_connected_components(graph)

    for idx, cc in enumerate(cc_list):
        cc_peptides = set()
        logger.debug(f"Connected component {idx}: {cc}")        
        for node, attributes in cc.nodes(data=True):
            nx.set_node_attributes(cc, {node: {'component_id': idx}})
            # TODO: adding segmement (locus) info as helpful format e.g. CM026586.1:734200-734540 

            peptides = attributes.get('peptides', None)
            if peptides is not None:
                cc_peptides.update(peptides)
        if is_only_one(cc_peptides):
            peptide_of_cc = get_only_one_peptide(cc_peptides, idx, input_peptides)            
            cc_with_only_one_peptide_list.append(peptide_of_cc)
            # update the graph with the epitope list (candidates)
            add_epitopes_on_cc(cc, [peptide_of_cc])
            logger.debug(f"Peptide only one for component {idx} - {peptide_of_cc.id}: {peptide_of_cc.seq}")
            
        else:  
            # logger.debug(f"Peptides on cc: {idx} : {cc_peptides}")
            peptides_on_cc_list = get_peptides(cc_peptides, idx, input_peptides)
            # add_epitopes_on_cc(cc, peptides_on_cc_list) # the epitope was not predicted yet
            peptides_on_cc_list.sort(key = lambda x: (len(x.seq), x.id), reverse=True)
            logger.debug(f"Number of peptides on cc: {idx} : {len(peptides_on_cc_list)}")
            logger.debug(f"List of peptides: {idx} : {peptides_on_cc_list}")
            _write_output_file_(peptides_on_cc_list, f"{cc_out_dir}/cc-{idx}.fasta") 
        

    cc_with_only_one_peptide_list.sort(key = lambda x: (len(x.seq), x.id), reverse=True)
    _write_output_file_(cc_with_only_one_peptide_list, f"{output_path}/cc_with_only_one_peptide.fasta")
    output_file_name = _get_output_file_name(graph_file, output_dir)
    # The graph was updated with the component_id and epitope_candidates from the cc with peptides lonely. 
    pickle.dump(graph, open(f"{output_file_name}-graph-cc-id.pickle", 'wb'))

def add_epitopes_on_cc(cc, epitope_candidates:list):
    for node, attributes in cc.nodes(data=True):
        nx.set_node_attributes(cc, {node: {'epitope_candidates': epitope_candidates}})

def get_new_seq_record(seq:SeqRecord, new_id:str):
    return SeqRecord(Seq(seq.seq), id=new_id, name="", description="")

def get_only_one_peptide(peptide_ids, idx, input_peptides):
    peptide_id = next(iter(peptide_ids))
    seq = input_peptides[peptide_id]
    new_id = f"{idx}-{peptide_id}"
    peptide = get_new_seq_record(seq, new_id)
    return peptide


def get_peptides(peptides_ids, idx, input_peptides):
    peptides_on_cc = []
    for pep_id in peptides_ids:
        peptide = input_peptides[pep_id]
        peptide.id = f"{idx}-{peptide.id}"
        peptide.name = ""
        peptide.description = ""
        peptides_on_cc.append(peptide)
    return peptides_on_cc

def is_only_one(peptides):
    return len(peptides) == 1

def identifying_and_sort_connected_components(graph):
    connected_component_list_result = []
    connected_component_list_sorted = sorted(nx.connected_components(graph), key=len, reverse=True)

    for cc in connected_component_list_sorted:
        connected_component_list_result.append(graph.subgraph(cc))

    return connected_component_list_result

def _get_output_path(input, outdir):
    outdir_path = os.path.dirname(input)
    if outdir:
        outdir_path = outdir
    return outdir_path

def _get_output_file_name(input, outdir):
    outdir_path = _get_output_path(input, outdir)

    input_file_name =  os.path.basename(input)
    input_file_name_splited =  os.path.splitext(input_file_name)
    output_file_base_name = input_file_name_splited[0]
    return os.path.join(outdir_path, output_file_base_name)

def _read_fasta(faa_file):
    return SeqIO.index(faa_file, "fasta")

def _write_output_file_(sequences, output_file_name):
    with open(output_file_name, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

if __name__ == "__main__":
    main()
