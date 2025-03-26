#!/usr/bin/env python

import click
import networkx as nx
import matplotlib.pyplot as plt
import pickle

@click.command(help="Aim of this program is to load and print a graph from a file.")
@click.argument("graph_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("connected_component", "--cc", type=int, default=None, help="The connected component to print.")
@click.option('quiet', '--q', is_flag=True, help="Do not print the graph information.")
def main(graph_file: click.Path, connected_component: int, quiet: bool):

        # Load the graph from the pickle file
    with open(graph_file, 'rb') as f:
        graph = pickle.load(f)
    
    # Print basic information about the graph
    print(f"Number of nodes: {graph.number_of_nodes()}")
    print(f"Number of edges: {graph.number_of_edges()}")
    
  
    cc_list = identifying_and_sort_connected_components(graph)
    
    for idx, cc in enumerate(cc_list):
        if connected_component is not None and idx == connected_component:
            if not quiet:
                print_cc_details(cc, idx)
            draw_graph(cc)
        elif connected_component is None and not quiet:
            print_cc_details(cc, idx)
    

    if connected_component is None:
        draw_graph(graph)

def identifying_and_sort_connected_components(graph):
    connected_component_list_result = []
    connected_component_list_sorted = sorted(nx.connected_components(graph), key=len, reverse=True)

    for cc in connected_component_list_sorted:
        connected_component_list_result.append(graph.subgraph(cc))

    return connected_component_list_result


def print_cc_details(cc, idx):
    print(f"Connected component {idx}: {cc}")        
    for node, attributes in cc.nodes(data=True):
        print(f"Node: {node}")
        for attr, value in attributes.items():
            print(f"  {attr}: {value}")    

def draw_graph(graph):
    plt.figure(figsize=(10, 10))
    pos = nx.spring_layout(graph)
    labels = {node: f"{node}\n{attributes.get('component_id', '')}" for node, attributes in graph.nodes(data=True)}
    nx.draw(graph, pos, labels=labels, node_size=500, font_size=10, font_weight="bold")
    plt.title("Graph Visualization")
    plt.show()


if __name__ == "__main__":
    main()
