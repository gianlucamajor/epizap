#!/usr/bin/env python

import logging.handlers
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
import logging


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

logger = logging.getLogger("SGC")
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

@click.command(help="Aim of this program is to create a Graph in witch the Nodes are the segments and the Edges indicate that there is a percentage (--threshold) of reads shared between two segments.")
@click.argument("segment_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--outdir", "-o", type=click.Path(exists=True, file_okay=False, dir_okay=True),  help="The dir path where the output file will be created.")
@click.option("--threshold", "-t", type=float, default=1, help="The min value [between 0.00 and 1] to create an EDGE between two NODES.")
@click.option("--processors", "-p", type=int, default=2)
def main(segment_file:click.Path, outdir:click.Path, threshold:float, processors:int):
   
    
    start_time = time.time()
    _setup_csv_field_size_limit()
    segments = []

    logger.info(f"threshold value setup: {threshold}")
    logger.info("starting read CSV file")

    
    line = get_next_line(segment_file)
    t = Pool(processes=processors)

    for l in line:
        for seg in  t.map(process_line, (l,)):
            segments.append(seg)
    t.close()
    t.join()

            
    logger.info("Ending of Read CSV file")
    logger.info("--- %s seconds ---" % (time.time() - start_time))
    logger.info("Starting sorting segments")
    logger.info("Starting sorting segments")

    segments.sort(key=Segment.get_complete_name, reverse=False)

    logger.info("Ending of sorting segments")
    logger.info("--- %s seconds ---" % (time.time() - start_time))

    logger.info("Starting graph handler")
    s_graph = nx.Graph()

    logger.info("adding nodes")
    for s in segments:
        s_graph.add_node(s.get_complete_name(), reads=s.get_set_of_reads(), mappings_number=s.get_mappings_number(), peptides=s.get_distinct_peptide_ids())

    logger.info("--- %s seconds ---" % (time.time() - start_time))
    logger.info("search and adding edges")
    segment_reads = {s.get_complete_name(): s.get_set_of_reads() for s in segments}
    for s_name, s_reads in segment_reads.items():
        for sc_name, sc_reads in segment_reads.items():
            if s_name != sc_name:
                set_shared = s_reads.intersection(sc_reads)
                p_shared = len(set_shared) / len(s_reads)
                if p_shared >= threshold:  # there is an edge
                    s_graph.add_edge(s_name, sc_name)
    
    logger.info("Graph: %s", s_graph)
    logger.info("--- %s seconds ---" % (time.time() - start_time))

    logger.info("saving graph on file")
    output_file_name = get_output_file_name(segment_file, outdir)
    
    pickle.dump(s_graph, open(f"{output_file_name}-graph.pickle", 'wb'))
    logger.info("--- %s seconds ---" % (time.time() - start_time))

    # logger.info("drawing graph")
    # nx.draw_networkx(s_graph)
    # plt.draw()
    # print("--- %s seconds ---" % (time.time() - start_time))
    # plt.show()


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
