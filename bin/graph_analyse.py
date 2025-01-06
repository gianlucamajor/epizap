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
@click.option("--segment", "-s", type=click.Path(exists=True, file_okay=True, dir_okay=False), default=None)
@click.option("--outdir", "-o", type=click.Path(exists=True, file_okay=False, dir_okay=True),  help="The dir path where the output file will be created.")
def main(graph_file:click.Path, epitopes_file:click.Path, segment:click.Path, outdir:click.Path):
    
    output_file = get_output_file_name(epitopes_file, outdir)
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

        # print(seg_inserts_id, ept_id)

        ept_seg = ept_by_segments.get(seg_inserts_id)
        if not ept_seg:
            epts = []
            
            ept_element = {}
            ept_element['id'] = ept_id
            ept_element['seq'] = str(seq)
            epts.append(ept_element)
            ept_by_segments[seg_inserts_id] = epts

        else: 
            ept_element = {}
            ept_element['id'] = ept_id
            ept_element['seq'] = str(seq)
            ept_seg.append(ept_element)
        
    # print(ept_by_segments.get("WNWZ01000402.1-22632-22735-2"))
    # print(ept_by_segments.get("CM026620.1-161873-162157-827"))
    # print(ept_by_segments.get("WNWZ01000039.1-97394-97488-252"))
    return ept_by_segments
        
        
    

def get_graph(graph_file, epitopes, output_file):
    G = pickle.load(open(graph_file, 'rb'))
    print("Graph:",G)

    
    
    # print(nx.clustering(G))
    # print(G.nodes())

    # for n in G.nodes():
    #     cc = nx.node_connected_component(G,n)
    #     print(n, len(cc))

    # cc = nx.node_connected_component(G, "CM026603.1-94787-94954-11")
    print(nx.number_connected_components(G))
    
    all_cc = nx.connected_components(G)

    # largest = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    # print(len(largest))

    S = []

    all_cc_sorted = sorted(nx.connected_components(G), key=len, reverse=False)
    # # # [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]

    for cc_s in all_cc_sorted:
        if len(cc_s) > 0:
            S.append(G.subgraph(cc_s))
        # print(len(cc_s))
    final_ept_list = []
    final_ept_set = set()
    for idx, sub_graph in enumerate(S):
        # each graph is a coonected componet
        # 2200
        # if(idx == 7000):
        # if "CM026617.1-197083-197351-497-27" in graph.nodes(): 
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
            
        # print(idx, epts)
        # if "CM026590.1-343373-343722-323" in graph.nodes(): 

        
        if idx == 5973:
            print(idx, epts_seq_set_of_graph, sub_graph.nodes())
            nx.draw_spring(sub_graph, with_labels=True)
            plt.show()
        print(idx, len(sub_graph.nodes()), len(sub_graph.edges()), len(epts_seq_set_of_graph))  
        unique_reads = set()
        for n in sub_graph.nodes(data=True):
            unique_reads.update(n[1].get('reads', []))
        print(f"Total unique reads in subgraph: {len(unique_reads)}")
        # for n in sub_graph.nodes(data=True):
        #     print(f"Node: {n[0]}, Reads: {len(n[1].get('reads'))}")

    
    total_of_epitopes = 0
    epitopes_sequence_final = []
    for idx_cc, conex_component_ept_list in enumerate(final_ept_list):
        for idx_segment, ept in enumerate(conex_component_ept_list):
            epitopes_sequence_final.append(_create_seq_record(ept, idx_cc,idx_segment))
            # print(idx_cc, idx_segment, ept)
        total_of_epitopes += len(conex_component_ept_list)
    
    #/home/gianluca/workspace/epizap/results_21_11_2024/epitopes/epitopes-final-list.fasta
    _write_output_file_(epitopes_sequence_final, output_file)

    # print(total_of_epitopes)
    print(total_of_epitopes)
    # print(len(final_ept_set))
    # print(len(final_ept_list))
       
   
    
    # nx.draw_spring(S[0], with_labels=True)
    # # plt.draw()
    # plt.show()

def get_output_file_name(input, outdir):
    outdir_path = os.path.dirname(input)
    if outdir:
        outdir_path = outdir

    # input_file_name =  os.path.basename(input)
    # input_file_name_splited =  os.path.splitext(input_file_name)
    # output_file_base_name = input_file_name_splited[0]
    output_file_base_name = "epitopes-cc-graph.fasta"
    return os.path.join(outdir_path, output_file_base_name) 

def _rename_ids_and_names(file, output):
    sequences = []
    list_of_sequence = _read_fasta(file)
    for seq_id in list_of_sequence:
        try:
            r = list_of_sequence[seq_id]
            id = r.id
            name = r.name
            r.id = id.replace("control_and_chagasic_patients_", "").replace("_", "-") 
            r.name = name.replace("control_and_chagasic_patients_", "").replace("_", "-") 
            r.description = ""
            sequences.append(r)
        except Exception as e:
            print("Ops", e, "occurred")
    _write_output_file_(sequences, output)
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

def _create_seq_record(sequence, name, idx):
    name_with_idx = f"{name}-{idx}"
    return SeqRecord(Seq(sequence), id=name_with_idx, name=name_with_idx, description="")

if __name__ == "__main__":
    main()
