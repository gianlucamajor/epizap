#!/usr/bin/env python

from model.Segment import Segment
import os 

import click
import networkx as nx
import matplotlib.pyplot as plt
import pickle

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

@click.command(help="Aim of this program is ...")
@click.argument("graph_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument("epitopes_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--outdir", "-o", type=click.Path(exists=True, file_okay=False, dir_okay=True),  help="The dir path where the output file will be created.")
def main(graph_file:click.Path, epitopes_file:click.Path, outdir:click.Path):
    
    output_file = _get_output_file_name(epitopes_file, outdir)
    epitopes = get_epitopes(epitopes_file)
    get_graph(graph_file, epitopes, output_file)



def get_epitopes(ept_file):
    epitopes_dict = _get_dic_from_fasta(ept_file)
    
    # >CM026586.1-734200-734540-29532-677-0
    # >control_and_chagasic_patients_CM026583.1-1048335-1048450-3-1-114361
    
    ept_by_segments = {}
    for k in epitopes_dict:
        complete_segment_id  = epitopes_dict[k].id
        seq = epitopes_dict[k].seq
        seg_inserts_and_pep_id = complete_segment_id.rsplit("-",1)[0] 
        ept_id = complete_segment_id.rsplit("-",1)[1]
        seg_inserts_id = seg_inserts_and_pep_id.rsplit("-",1)[0]

        # print(complete_segment_id, seg_inserts_and_pep_id, seg_inserts_id, ept_id)

        ept_seg = ept_by_segments.get(seg_inserts_and_pep_id)
        if not ept_seg:
            epts = []
            
            ept_element = {}
            ept_element['id'] = ept_id
            ept_element['seq'] = str(seq)
            epts.append(ept_element)
            ept_by_segments[seg_inserts_and_pep_id] = epts

        else: 
            ept_element = {}
            ept_element['id'] = ept_id
            ept_element['seq'] = str(seq)
            ept_seg.append(ept_element)
        
    return ept_by_segments
        
        
    

def get_graph(graph_file, epitopes, output_file):
    G = pickle.load(open(graph_file, 'rb'))
    # print("Graph:",G)
    # print(nx.number_connected_components(G))

    print(F"CC idx \t Nodes \t Edges \t Epitopes \t Unique Reads")
    
    all_cc = nx.connected_components(G)

    S = []

    all_cc_sorted = sorted(nx.connected_components(G), key=len, reverse=False)

    for cc_s in all_cc_sorted:
        if len(cc_s) > 0:
            S.append(G.subgraph(cc_s))
            
    final_ept_list = []
    final_ept_set = set()
    total_of_unique_reads_on_cc = 0
    for idx, sub_graph in enumerate(S):
        epts = []
        epts_seq_set_of_graph = set()
        for n in sub_graph.nodes():
            epts_from_segment = epitopes.get(n)
            epts.append(epts_from_segment)

            if epts_from_segment:
                for e in epts_from_segment:
                    epts_seq_set_of_graph.add(e['seq'])
                    final_ept_set.add(e['seq'])
    
        final_ept_list.append(epts_seq_set_of_graph)
        
        ## START - print the number of nodes, edges and epitopes in the connected component
        """
        if idx == 6815:

            unique_reads = set()
            for n in sub_graph.nodes(data=True):
                print("len:", len(n[1].get('reads', [])))
                unique_reads.update(n[1].get('reads', []))

            # To create a list of links to IGV  browser of segments and its insertions
            igv_links = []
            for n in sub_graph.nodes():
                igv_links.append(_create_igv_link(n))
            igv_links_str = ", ".join(igv_links)

            print(idx, epts_seq_set_of_graph, sub_graph.nodes())
            # print(idx, len(sub_graph.nodes()), len(sub_graph.edges()), len(epts_seq_set_of_graph))  
            print(F"{idx} \t {len(sub_graph.nodes())} \t {len(sub_graph.edges())} \t {len(epts_seq_set_of_graph)} \t {len(unique_reads)} \t {igv_links} \t {sub_graph.nodes()}")  
            nx.draw_spring(sub_graph, with_labels=True)
            plt.show()
        """
        ## END - print the number of nodes, edges and epitopes in the connected component

        
        unique_reads = set()
        for n in sub_graph.nodes(data=True):
            unique_reads.update(n[1].get('reads', []))
        # CC index, number of nodes, number of edges, number of epitopes, number of unique reads
        
        # To create a list of links to IGV  browser of segments and its insertions and MView
        igv_links = []
        mview_links = []
        for n in sub_graph.nodes():
            igv_links.append(_create_igv_link(n))
            mview_links.append(_create_mview_link(n))
        igv_links_str = ", ".join(igv_links)
        mview_links_str = ", ".join(mview_links)


        print(F"{idx} \t {len(sub_graph.nodes())} \t {len(sub_graph.edges())} \t {len(epts_seq_set_of_graph)} \t {len(unique_reads)} \t {', '.join(epts_seq_set_of_graph)} \t {igv_links_str} \t {mview_links_str}")
        
        
        # print(f"Total unique reads in subgraph: {len(unique_reads)}")
        total_of_unique_reads_on_cc += len(unique_reads) 
        # for n in sub_graph.nodes(data=True):
        #     print(f"Node: {n[0]}, Reads: {len(n[1].get('reads'))}")

    
    total_of_epitopes = 0
    epitopes_sequence_final = []
    for idx_cc, conex_component_ept_list in enumerate(final_ept_list):
        for idx_segment, ept in enumerate(conex_component_ept_list):
            epitopes_sequence_final.append(_create_seq_record(ept, idx_cc,idx_segment))
        total_of_epitopes += len(conex_component_ept_list)
    
    _write_output_file_(epitopes_sequence_final, output_file)

    # print(total_of_epitopes)
    print(total_of_epitopes)
    print(total_of_unique_reads_on_cc)
    # print(len(final_ept_set))
    # print(len(final_ept_list))

def _create_mview_link(segment):
    local_host = "http://localhost:8080"
    lbi_host = "http://projetos.lbi.iq.usp.br/trypanosoma/epizap"
    mview_url = "mview/"
   
    base_url = f"{lbi_host}/{mview_url}"

    parts = segment.split('-')
    if len(parts) < 5:
        raise ValueError("Invalid segment ID format")
    
    number_of_peptides = int(parts[4])
    if number_of_peptides <= 1:
        return "There is no MSA"
    
    msa_id = f"{parts[0]}-{parts[1]}-{parts[2]}-{parts[3]}-{number_of_peptides}"
    return f"{base_url}{msa_id}.html"


def _create_igv_link(segment_id):
    """
        Expected is received something like this:
        CM026600.1-29631-29710-758
        And return a link like this:
        http://localhost:8080/igv-webapp/?locus=CM026600.1:29631-29710
    """

    local_host = "http://localhost:8080"
    lbi_host = "http://projetos.lbi.iq.usp.br/trypanosoma/epizap"
    igv_webapp = "igv-webapp/?locus="

    base_url = f"{lbi_host}/{igv_webapp}"
    
    
    parts = segment_id.split('-')
    if len(parts) < 3:
        raise ValueError("Invalid segment ID format")
    locus = f"{parts[0]}:{parts[1]}-{parts[2]}"
    return f"{base_url}{locus}"
 
def _get_output_file_name(input, outdir):
    outdir_path = os.path.dirname(input)
    if outdir:
        outdir_path = outdir

    output_file_base_name = "epitopes-cc-graph.fasta"
    return os.path.join(outdir_path, output_file_base_name) 


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

def _create_seq_record(sequence, name, idx):
    name_with_idx = f"{name}-{idx}"
    return SeqRecord(Seq(sequence), id=name_with_idx, name=name_with_idx, description="")

if __name__ == "__main__":
    main()
