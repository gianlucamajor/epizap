#!/usr/bin/env python

import click
import logging
import pickle
import networkx as nx
import csv
import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# from model.Segment import Segment
from AnnotaionHandler import AnnotationHandler

logger = logging.getLogger("GU")
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

CC_ID_POS = 0
EPITOPES_POS = 2

@click.command(help="Aim of this program is ...")
@click.argument("graph_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument("epitopes_report_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--features_file", "-f", type=click.Path(exists=True, file_okay=True, dir_okay=False), default=None,  help="The file with features annotated on the segment.")
@click.option("--outdir", "-o", type=click.Path(exists=True, file_okay=False, dir_okay=True),  help="The dir path where the output file will be created.") 
def main(graph_file:click.Path, epitopes_report_file:click.Path, features_file:click.Path, outdir:click.Path):

    # cc_id	number_of_epitopes	epitopes	consensus	msa_avg_identity
    epitope_candidates_reported_dic = {}
    with open(epitopes_report_file, "r") as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        next(reader)  # Skip the header line
        for row in reader:
            epitope_candidates_reported_dic[int(row[CC_ID_POS].strip())] = row
    
    if features_file:
        aFeaturesHandler = AnnotationHandler(features_file)
        


    with open(graph_file, 'rb') as f:
        graph = pickle.load(f)

    # Print basic information about the graph
    logger.info(f"Number of nodes: {graph.number_of_nodes()}")
    logger.info(f"Number of edges: {graph.number_of_edges()}")

    cc_sorted_list = get_sorted_connected_components(graph)
    
    for cc in cc_sorted_list:
        logger.debug(f"Connected component: {cc}")
        for node, attributes in cc.nodes(data=True):
            cc_ids = set()
            logger.debug(f"Node: {node}")
            if features_file:
                feature = _get_features_annotation(node, aFeaturesHandler)
                add_feature_on_cc(cc, node, feature)

            for attr, value in attributes.items():                
                if attr == 'component_id':
                    cc_ids.add(value)
                
            if len(cc_ids) > 1:
                raise Exception(f"CC {cc} has multiple component_ids: {cc_ids}")
        cc_id = cc_ids.pop()
        logger.debug(f"component_id: {cc_id}")
        epitope_candidate_report = epitope_candidates_reported_dic.get(cc_id, None)
        ## The epitopes that don't require MSA, such cc where there is only one peptide, already should be in the graph.
        if epitope_candidate_report is not None:
            add_epitopes_on_cc(cc, eval(epitope_candidate_report[EPITOPES_POS]))
            # seq_record_list = eval(epitope_candidates_list)
            # print(seq_record_list)
    
    output_file_name = _get_output_file_name(graph_file, outdir)
    pickle.dump(graph, open(f"{output_file_name}-msa-epitopes.pickle", 'wb'))

def _get_features_annotation(segment_id, annHandler):
    
    parts = segment_id.split('-')
    if len(parts) < 3:
        raise ValueError("Invalid segment ID format")
    
    scaffold = parts[0]
    start = int(parts[1])
    end = int(parts[2])

    return annHandler.get_feature(scaffold, start, end, "CDS")

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


def add_epitopes_on_cc(cc, epitope_candidates:list):
    for node, attributes in cc.nodes(data=True):
        nx.set_node_attributes(cc, {node: {'epitope_candidates': epitope_candidates}})

def add_feature_on_cc(cc, node, feature):
    nx.set_node_attributes(cc, {node: {'feature': feature}})

def get_sorted_connected_components(graph):
    connected_component_list_result = []
    connected_component_list_sorted = sorted(nx.connected_components(graph), key=len, reverse=True)

    for cc in connected_component_list_sorted:
        connected_component_list_result.append(graph.subgraph(cc))

    return connected_component_list_result

if __name__ == "__main__":
    main()
