#!/usr/bin/env python

import sys
import csv
import os
import time

from model.Segment import Segment
from multiprocessing.pool import ThreadPool as Pool

import click
import networkx as nx
import matplotlib.pyplot as plt
import pickle


TAB_DELIM="\t"
SCF_POS=0
START_POS=1
END_POS=2
TOTAL_READS_MAPPED_POS=3
MIN_MAPQ=4
MAX_MAPQ=5
AVG_MAPQ=6
MEDIAN_MAPQ=7
LIST_OF_READS=8

@click.command(help="Aim of this program is compute adjacency matrix from a list of segments...")
@click.argument("segment_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--outdir", "-o", type=click.Path(exists=True, file_okay=False, dir_okay=True),  help="The dir path where the output file will be created.")
@click.option("--threshold", "-t", type=float, default=1, help="The min value [between 0.00 and 1] to create an EDGE between two NODES.")
@click.option("--processors", "-p", type=int, default=2)
def main(segment_file:click.Path, outdir:click.Path, threshold:float, processors:int):
    start_time = time.time()
    _setup_csv_field_size_limit()
    segments = []
    
    print("staring read CSV file")

    
    line = get_next_line(segment_file)
    t = Pool(processes=processors)

    for l in line:
        for seg in  t.map(process_line, (l,)):
            segments.append(seg)
    t.close()
    t.join()

            
    print("Ending of Read CSV file")
    print("--- %s seconds ---" % (time.time() - start_time))

    print("Starting sorting segments")

    segments.sort(key=Segment.get_name, reverse=False)
    segments_copy = segments.copy()

    print("Ending of sorting segments")
    print("--- %s seconds ---" % (time.time() - start_time))

    print("Starting graph handler")
    s_graph = nx.Graph()

    print("adding nodes")
    for s in segments:
        s_graph.add_node(s.get_name())

    print("--- %s seconds ---" % (time.time() - start_time))
    print("search and adding edges")
    for idx, s in enumerate(segments):
        for sc in segments_copy:
            set_shared = s.get_set_of_reads().intersection(sc.get_set_of_reads())
            p_shared = len(set_shared) / len(s.get_set_of_reads())
            # print(p_shared)
            if s.get_name() != sc.get_name() and p_shared >= threshold: # there is a edge
                s_graph.add_edge(s.get_name(), sc.get_name())
    
    print(":", s_graph)
    print("--- %s seconds ---" % (time.time() - start_time))

    print("saving graph on file")
    # output_file_name = get_output_file_name(segment_file)
    output_file_name = get_output_file_name(segment_file, outdir)
    
    pickle.dump(s_graph, open(f"{output_file_name}-graph.pickle", 'wb'))
    print("--- %s seconds ---" % (time.time() - start_time))

    print("drawing graph")
    nx.draw_networkx(s_graph)
    plt.draw()
    print("--- %s seconds ---" % (time.time() - start_time))
    plt.show()

def get_next_line(segment_file):
    with open(segment_file) as file:
        file_reader = csv.reader(file, delimiter=TAB_DELIM)
        for line in file_reader:
            yield line

def process_line(line):
    return Segment(line[SCF_POS], line[START_POS], line[END_POS], line[TOTAL_READS_MAPPED_POS] ,line[MIN_MAPQ], line[MAX_MAPQ], line[AVG_MAPQ], line[MEDIAN_MAPQ], line[LIST_OF_READS].split(";") )

def _setup_csv_field_size_limit():
    # increasing csv file size limit
    csv.field_size_limit(int(sys.maxsize/1000))

def get_output_file_name(input, outdir):
    outdir_path = os.path.dirname(input)
    if outdir:
        outdir_path = outdir

    input_file_name =  os.path.basename(input)
    input_file_name_splited =  os.path.splitext(input_file_name)
    output_file_base_name = input_file_name_splited[0]
    return os.path.join(outdir_path, output_file_base_name) 
    

if __name__ == "__main__":
    main()
