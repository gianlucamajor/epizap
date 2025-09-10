
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
from BlastResults import BlastResults
from IEDBEpitopeTableHandler import IEDBEpitopeTableHandler

logger = logging.getLogger("ER")
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


@click.command(help="Aim of this program is to generate a report of the epitopes from the graph.")
@click.argument("graph_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--single-reads/--no-single-reads", "single_reads", default=False, help="sequence from CC composed by single reads will not be included in the fasta file reported by default.")
@click.option("--iedb", "-iedb", is_flag=True, help="Include IEDB information in the report.")
@click.option("--iedb-epitopes", "-iedb-path", "ept_file_path", type=click.Path(exists=True, file_okay=True, dir_okay=False),  help="The path to the IEDB epitope table CSV file. Required if --iedb is set.")
@click.option("--iedb-blast-hits", "-iedb-hits-path", "iedb_blast_file_path", type=click.Path(exists=True, file_okay=True, dir_okay=False),  help="The path to the BLAST results file of the epitopes against the IEDB epitope database. Required if --iedb is set.")
@click.option("--outdir", "-o", type=click.Path(exists=True, file_okay=False, dir_okay=True),  help="The dir path where the output file will be created.")
@click.option('-v', '--verbose', is_flag=True)
def main(graph_file:click.Path, 
         outdir:click.Path, 
         single_reads:bool, 
         iedb:bool, 
         ept_file_path:click.Path, 
         iedb_blast_file_path:click.Path, 
         verbose:bool):

    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose mode is on.")

    with open(graph_file, 'rb') as f:
        graph = pickle.load(f)

    # Print basic information about the graph
    logger.info(f"Number of nodes: {graph.number_of_nodes()}")
    logger.info(f"Number of edges: {graph.number_of_edges()}")
    logger.info(f"IEDB Info: {iedb}")
    logger.info(f"Single Reads Allowed: {single_reads}")
    logger.info(f"IEDB Epitopes reported: {iedb}")
    
    validate_iedb_options(iedb, ept_file_path, iedb_blast_file_path)
    if iedb:
         blastIEDBHandler = BlastResults(iedb_blast_file_path)
         min_length = 8
         min_identity = 100
         no_gaps = True
         qseqid = None
         logger.info(f"IEDB Parameters: min_length: {min_length}, min_identity: {min_identity}, no_gaps: {no_gaps}, qseqid: {qseqid}")
         iedbEpitopesBestHits = blastIEDBHandler.best_hits_by_sseqid(min_length, min_identity, no_gaps, qseqid)
         IEDBTableHandler = IEDBEpitopeTableHandler(ept_file_path)

    
    logger.debug(f"CC_idx\tNof_Nodes\tNof_Edges\tNof_Epitopes\tNof_Unique_peptides\tNof_Unique_Reads\tEpitope_candidates\tMSA_links\tIGV_links\tFeatures")

    cc_sorted_list = get_sorted_connected_components(graph)
    epitope_candidates_graph = []
    report_json_epitopes_list = []
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
            for attr, value in attributes.items():
                if attr == 'component_id':
                    cc_ids.add(value)
                if attr == 'reads':
                    reads_cc.update(value)
                if attr == 'peptides':
                    pepiteds_cc.update(value)
                if attr == 'epitope_candidates' and len(epitope_candidates_cc) == 0: # all the nodes should has the same epitope candidates
                    epitope_candidates_cc = value
                if attr == 'feature':
                    features_cc.extend(parse_feature_string(value))

            if len(cc_ids) > 1:
                raise Exception(f"CC {cc} has multiple component_ids: {cc_ids}")
            
        cc_id = cc_ids.pop()

        msa_links_cc_= _create_mview_link(cc_id, pepiteds_cc)
        msa_page_name = _create_mview_name(cc_id, pepiteds_cc)

        if single_reads:
            epitope_candidates_graph.extend(epitope_candidates_cc)
        else:
            if len(reads_cc) >= 2: 
                epitope_candidates_graph.extend(epitope_candidates_cc)
            
        
        
        ept_candidantes_seq = []
        for e in epitope_candidates_cc:
            ept_candidantes_seq.append(str(e.seq))
            logger.debug(f"{e.id}\t{cc.number_of_nodes()}\t{cc.number_of_edges()}\t{len(epitope_candidates_cc)}\t{len(pepiteds_cc)}\t{len(reads_cc)}\t{str(e.seq)}\t{msa_links_cc_}\t{igv_links}\t{features_cc}")
            
            iedbEpitopesInfo = []
            if iedb:
               for iedbHits in iedbEpitopesBestHits:
                    epitope_pos = 0
                    iedbEpitope_pos = 1
                    qstart_pos = 6
                    qend_pos = 7 
                    sstart_pos = 8 
                    send_pos = 9
                

                    if(iedbHits[epitope_pos] == e.id):
                        iedbEpitopeInfo = IEDBTableHandler.get_by_epitope_id(iedbHits[iedbEpitope_pos])
                        logger.info(f"IEDB Best hit for epitope {iedbHits}")
                        logger.info(f"{iedbEpitopeInfo['Epitope_id']}\t{iedbEpitopeInfo['epitope_sequence']}\t{iedbEpitopeInfo['Epitope - Source Molecule']}\t{iedbEpitopeInfo['Epitope - Source Molecule IRI']}")
                        iedbEpitopesInfo.append({"IEDB_id": iedbEpitopeInfo['Epitope_id'],
                                               "sequence": iedbEpitopeInfo['epitope_sequence'],
                                               "qstart": iedbHits[qstart_pos],
                                               "qend": iedbHits[qend_pos],
                                               "sstart": iedbHits[sstart_pos],
                                               "send": iedbHits[send_pos],
                                               "source_molecule": iedbEpitopeInfo['Epitope - Source Molecule'],
                                               "source_molecule_IRI": iedbEpitopeInfo['Epitope - Source Molecule IRI']})

            # Create a JSON object with the same information
            json_data = {
                "ID": e.id,
                "Number of Genomic Regions": cc.number_of_nodes(),
                "Number of Peptides": len(pepiteds_cc),  # unique peptides
                "Number of Inserts": len(reads_cc),  # unique reads
                "Epitope": str(e.seq),
                "MSA": msa_page_name,
                "Genomic Region Locus": genomic_region_locus,
                "Features": {
                    "GenomicRegionsAnnotation": features_cc,
                    "TCruziIEDB": iedbEpitopesInfo,
                }
            }
            report_json_epitopes_list.append(json_data)

    
    output_file_name = _get_output_file_name(graph_file, outdir)
    _write_output_fasta_file(epitope_candidates_graph, f"{output_file_name}.fasta")
    
    # Write the JSON object to a file
    with open(f"{output_file_name}.json", "w") as json_file:
        json.dump(report_json_epitopes_list, json_file, indent=4)
        
def validate_iedb_options(iedb, ept_file_path, iedb_blast_file_path):
    if iedb:
        if not ept_file_path or not iedb_blast_file_path:
            logger.error(
                "When --iedb is set, you must also provide both --iedb-epitopes/-epitopes and --iedb-blast-hits/-epitopes-hits options."
            )
            raise click.UsageError(
                "Missing required options: --iedb-epitopes/-epitopes and/or --iedb-blast-hits/-epitopes-hits when --iedb is used."
            )
    
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

if __name__ == "__main__":
    main()