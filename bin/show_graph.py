#!/usr/bin/env python

import click
import networkx as nx
import matplotlib.pyplot as plt
import pickle

@click.command(help="Aim of this program is to load and print a graph from a file.")
@click.argument("graph_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("connected_component", "-cc", type=int, default=None, help="The connected component to print.")
@click.option('verbose', '-v', is_flag=True, help="Print detailed information about the connected components.")
def main(graph_file: click.Path, connected_component: int, verbose: bool):

        # Load the graph from the pickle file
    with open(graph_file, 'rb') as f:
        graph = pickle.load(f)
    
    # Print basic information about the graph
    print(f"Number of nodes: {graph.number_of_nodes()}")
    print(f"Number of edges: {graph.number_of_edges()}")
    
  
    cc_list = identifying_and_sort_connected_components(graph)
    
    for idx, cc in enumerate(cc_list):
        if connected_component is not None and idx == connected_component:
            print_cc(cc, idx, verbose)
            # draw_graph(cc)
        elif connected_component is None:
            print_cc(cc, idx, verbose)
    

    if connected_component is None:
        draw_graph(graph)

def identifying_and_sort_connected_components(graph):
    connected_component_list_result = []
    connected_component_list_sorted = sorted(nx.connected_components(graph), key=len, reverse=True)

    for cc in connected_component_list_sorted:
        connected_component_list_result.append(graph.subgraph(cc))

    return connected_component_list_result

def print_cc(cc, idx, verbose):
    if verbose:
        print_cc_details(cc, idx)
    else:
        print_cc_summary_in_line(cc, idx)
        # print_cc_summary(cc, idx)

def print_cc_summary(cc, idx):
    print(f"Connected component: {idx} - {cc}")        
    for node, attributes in cc.nodes(data=True):
        print(f"Node: {node}")
        for attr, value in attributes.items():
            if attr == 'peptides' or attr == 'reads':
                value = len(value)
            
            print(f"  {attr}: {value} ")

def print_cc_summary_in_line(cc, idx):
    for node, attributes in cc.nodes(data=True):
        line_array = []
        line_array.append(node)
        # print(f"{node}")        
        for attr, value in attributes.items():
            if attr == 'peptides' or attr == 'reads':
                value = len(value)
            line_array.append(value)
            # print(f"{attr}:\t{value}", end=" ")
        print("\t".join(str(item) for item in line_array))

def print_cc_details(cc, idx):
    print(f"Connected component {idx}: {cc}")        
    for node, attributes in cc.nodes(data=True):
        print(f"Node: {node}")
        for attr, value in attributes.items():
            print(f"  {attr}: {value}")    

def draw_graph(graph):
    plt.figure(figsize=(8, 8))
    pos = nx.spring_layout(graph, seed=42)  # Fixed seed for consistent layout
    labels = {node: f"{node}\n{attributes.get('component_id', '')}" for node, attributes in graph.nodes(data=True)}
    
    # Draw nodes with color based on attributes (if available)
    node_colors = [attributes.get('color', 'skyblue') for _, attributes in graph.nodes(data=True)]
    nx.draw_networkx_nodes(graph, pos, node_color=node_colors, node_size=1000, alpha=0.9)
    
    # Draw edges with transparency
    nx.draw_networkx_edges(graph, pos, alpha=0.5, edge_color="gray")
    
    # Draw labels
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=10, font_weight="bold")
    
    plt.title("Graph Visualization", fontsize=16)
    plt.axis("off")  # Turn off the axis
    plt.tight_layout(pad=0.5)  # Adjust layout to avoid clipping and increase padding
    plt.show()


if __name__ == "__main__":
    main()
