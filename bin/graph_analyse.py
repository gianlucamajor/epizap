#!/usr/bin/env python

from model.Segment import Segment

import click
import networkx as nx
import matplotlib.pyplot as plt
import pickle

@click.command(help="Aim of this program is ...")
def main():
    G = pickle.load(open('/home/gianluca/workspace/epizap/results_05_11_2024/segments/control_and_chagasic_patients-segments-on-same-strand-all.pickle', 'rb'))

    # G = pickle.load(open('/home/gianluca/workspace/epizap/results/segments/asympto_b_mapped-segment.pickle', 'rb'))
    print(G)
    
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
        if len(cc_s) > 1:
            S.append(G.subgraph(cc_s))
        print(len(cc_s))

    for idx, graph in enumerate(S):
        if(idx == 2101):
        # if "CM026592.1-138229-138326-6" in graph.nodes():
            print(graph.nodes())
            nx.draw_spring(graph, with_labels=True)
            plt.show()
        print(idx, len(graph.nodes()))
   
    
    # nx.draw_spring(S[0], with_labels=True)
    # # plt.draw()
    # plt.show()

if __name__ == "__main__":
    main()
