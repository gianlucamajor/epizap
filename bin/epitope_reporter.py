
#!/usr/bin/env python

import json
import logging
import click
import networkx as nx
import pickle
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger("ER")
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


@click.command(help="Aim of this program is to generate a report of the epitopes from the graph.")
@click.argument("graph_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
# @click.argument("epitopes_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--single_reads_not_allowed", "-no-single-reads", is_flag=True, help="sequence from CC composed by single reads will not be included in the fasta file reported.")
@click.option("--outdir", "-o", type=click.Path(exists=True, file_okay=False, dir_okay=True),  help="The dir path where the output file will be created.")
def main(graph_file:click.Path, outdir:click.Path, single_reads_not_allowed:bool):

    with open(graph_file, 'rb') as f:
        graph = pickle.load(f)

    # Print basic information about the graph
    logger.debug(f"Number of nodes: {graph.number_of_nodes()}")
    logger.debug(f"Number of edges: {graph.number_of_edges()}")
    
    print(f"CC_idx\tNof_Nodes\tNof_Edges\tNof_Epitopes\tNof_Unique_peptides\tNof_Unique_Reads\tEpitope_candidates\tMSA_links\tIGV_links\tFeatures")

    cc_sorted_list = get_sorted_connected_components(graph)
    epitope_candidates_graph = []
    report_json_epitopes_list = []
    report_epitopes_list = []
    for cc in cc_sorted_list:
        logger.debug(f"Connected component: {cc}")        
        reads_cc = set()
        pepiteds_cc = set()
        features_cc = []
        epitope_candidates_cc = []
        igv_links = []
        genomic_region_locus = []
        for node, attributes in cc.nodes(data=True):
            cc_ids = set()
            genomic_region_locus.append(_get_locus_from_segment(node))
            igv_links.append(_create_igv_link(node))
            # print(f"Node: {node}")
            for attr, value in attributes.items():
                # print(f"  {attr}: {value}")
                if attr == 'component_id':
                    cc_ids.add(value)
                if attr == 'reads':
                    reads_cc.update(value)
                if attr == 'peptides':
                    pepiteds_cc.update(value)
                if attr == 'epitope_candidates' and len(epitope_candidates_cc) == 0: # all the nodes should has the same epitope candidates
                    epitope_candidates_cc = value
                    # print(f"  {attr}: {value}")
                if attr == 'feature':
                    features_cc.extend(parse_feature_string(value))

            if len(cc_ids) > 1:
                raise Exception(f"CC {cc} has multiple component_ids: {cc_ids}")
            
        cc_id = cc_ids.pop()

        msa_links_cc_= _create_mview_link(cc_id, pepiteds_cc)
        msa_page_name = _create_mview_name(cc_id, pepiteds_cc)

        if single_reads_not_allowed:
            if len(reads_cc) >= 2: 
                epitope_candidates_graph.extend(epitope_candidates_cc)
        else:
            epitope_candidates_graph.extend(epitope_candidates_cc)
        
        
        ept_candidantes_seq = []
        for e in epitope_candidates_cc:
            ept_candidantes_seq.append(str(e.seq))
            # CC_idx	Nof_Nodes	Nof_Edges	Nof_Epitopes    Nof_Unique_peptides     Nof_Unique_Reads	Epitope_candidates     MSA_links    IGV_links   Features
            # print(f"{cc_id}\t{cc.number_of_nodes()}\t{cc.number_of_edges()}\t{len(epitope_candidates_cc)}\t{len(pepiteds_cc)}\t{len(reads_cc)}\t{ept_candidantes_seq}\t{msa_links_cc_}\t{igv_links}\t{features_cc}")
            # Print the information in tabular format
            print(f"{e.id}\t{cc.number_of_nodes()}\t{cc.number_of_edges()}\t{len(epitope_candidates_cc)}\t{len(pepiteds_cc)}\t{len(reads_cc)}\t{str(e.seq)}\t{msa_links_cc_}\t{igv_links}\t{features_cc}")

            # Create a JSON object with the same information
            json_data = {
                "ID": e.id,
                "Number of Genomic Regions": cc.number_of_nodes(),
                "Number of Peptides": len(pepiteds_cc), # unique peptides
                "Number of Inserts": len(reads_cc), #unique reads
                "Epitope": str(e.seq),
                "MSA": msa_page_name,
                "Genomic Region Locus": genomic_region_locus,
                "Features": features_cc
            }
            report_json_epitopes_list.append(json_data)

    
    output_file_name = _get_output_file_name(graph_file, outdir)
    _write_output_fasta_file(epitope_candidates_graph, f"{output_file_name}.fasta")
    
    # Write the JSON object to a file
    with open(f"{output_file_name}.json", "w") as json_file:
        json.dump(report_json_epitopes_list, json_file, indent=4)
        

def parse_feature_string(features_str):
    # Example input: ["WNWZ01000121.1-8261-8334 | KAF8303398.1 | tryptophanyl-tRNA synthetase, putative | 1.0", "WNWZ01000189.1-16368-16457 | pseudogene | tryptophanyl-tRNA synthetase, pseudonote | 0.8640776699029126"]
    features = []
    for feat in features_str: 
        if feat:
            parts = [p.strip() for p in feat.split('|')]
            if len(parts) != 4:
                raise ValueError(f"Invalid feature string format: {feat}")
            features.append({
        "genomic_region": parts[0],
        "type": parts[1],
        "description": parts[2],
        "coverage": float(parts[3])
            })
    
    return features

def get_sorted_connected_components(graph):
    connected_component_list_result = []
    connected_component_list_sorted = sorted(nx.connected_components(graph), key=len, reverse=True)

    for cc in connected_component_list_sorted:
        connected_component_list_result.append(graph.subgraph(cc))

    return connected_component_list_result

def _create_mview_name(cc_id, number_of_peptides):
    # Lonely peptides does have a MSA link
    mview_name = ""
    if len(number_of_peptides) > 1:
        mview_name = f"cc-{cc_id}.html"
    
    return mview_name

def _create_mview_link(cc_id, number_of_peptides):
    # Lonely peptides does have a MSA link
    mview_name = ""
    if len(number_of_peptides) > 1:
        local_host = "http://localhost:8080"
        lbi_host = "http://projetos.lbi.iq.usp.br/trypanosoma/epizap"
        mview_url = "mview/"
    
        base_url = f"{local_host}/{mview_url}"
        msa_id = f"cc-{cc_id}"
        mview_name = f"{base_url}{msa_id}.html"
    return mview_name

def _create_igv_link(node_name):
    """
        Expected is received something like this:
        CM026586.1-1665418-1665512-182-22
        And return a link like this:
        http://localhost:8080/igv-webapp/?locus=CM026586.1:1665418-1665512
    """

    local_host = "http://localhost:8080"
    lbi_host = "http://projetos.lbi.iq.usp.br/trypanosoma/epizap"
    igv_webapp = "igv-webapp/?locus="

    base_url = f"{local_host}/{igv_webapp}"
    locus = _get_locus_from_segment(node_name)
    return f"{base_url}{locus}"

def _get_locus_from_segment(segment_id):
    parts = segment_id.split('-')
    if len(parts) < 3:
        raise ValueError("Invalid segment ID format")
    return f"{parts[0]}:{parts[1]}-{parts[2]}"

def _write_output_fasta_file(sequences, output_file_name):
    with open(output_file_name, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

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

if __name__ == "__main__":
    main()