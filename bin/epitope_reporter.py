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
from BlastResults import BlastColumns
from IEDBEpitopeTableHandler import IEDBEpitopeTableHandler
from enum import Enum

class DefaultHistsParams(Enum):
    MIN_LENGTH = 8
    MIN_IDENTITY = 100
    NO_GAPS = True
    QSEQID = None

class ProteomeHistsParams(Enum):
    MIN_LENGTH = 8
    MIN_IDENTITY = 60
    NO_GAPS = False
    QSEQID = None
    
logger = logging.getLogger("ER")
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


@click.command(help="Aim of this program is to generate a report of the epitopes from the graph.")
@click.argument("graph_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--single-reads/--no-single-reads", "single_reads", default=False, help="sequence from CC composed by single reads will not be included in the fasta file reported by default.")
@click.option("--proteome-hits", "-proteome-hits", "proteome_hits_fp", type=click.Path(exists=True, file_okay=True, dir_okay=False),  help="The path to the BLAST results file of the epitopes against the T. cruzi proteome.")
@click.option("--iedb", "-iedb", is_flag=True, help="Include IEDB information in the report.")
@click.option("--iedb-tcruzi-epitopes-hits", "-tcruzi-epitopes-hits", "iedb_tcruzi_epitopes_hits_fp", type=click.Path(exists=True, file_okay=True, dir_okay=False),  help="The path to the BLAST results file of the epitopes against the IEDB epitope database. Required if --iedb is set.")
@click.option("--iedb-tcruzi-epitopes", "-tcruzi-epitopes", "iedb_tcruzi_epitopes_fp", type=click.Path(exists=True, file_okay=True, dir_okay=False),  help="The path to the IEDB T. cruzi epitope table CSV file. Required if --iedb is set.")
@click.option("--iedb-human-epitopes-hits", "-human-epitopes-hits", "iedb_human_epitopes_hits_fp", type=click.Path(exists=True, file_okay=True, dir_okay=False),  help="The path to the BLAST results file of the epitopes against the human IEDB epitope database. Optional.")
@click.option("--iedb-human-epitopes", "-human-epitopes", "iedb_human_epitopes_fp", type=click.Path(exists=True, file_okay=True, dir_okay=False),  help="The path to the IEDB human epitope table CSV file.")
@click.option("--inserts-group", "-inserts_group", "inserts_group_fp", type=click.Path(exists=True, file_okay=True, dir_okay=False),  help="The path to the inserts group file.")
@click.option("--outdir", "-o", type=click.Path(exists=True, file_okay=False, dir_okay=True),  help="The dir path where the output file will be created.")
@click.option('-v', '--verbose', is_flag=True)
def main(graph_file:click.Path, 
         outdir:click.Path, 
         single_reads:bool, 
         proteome_hits_fp:click.Path,
         iedb:bool, 
         iedb_tcruzi_epitopes_hits_fp:click.Path,
         iedb_tcruzi_epitopes_fp:click.Path, 
         iedb_human_epitopes_hits_fp:click.Path,
         iedb_human_epitopes_fp:click.Path,
         inserts_group_fp:click.Path,
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
    logger.info(f"IEDB  Epitopes report: {iedb}")
    logger.info(f"Inserts Group: {inserts_group_fp}")
    
    ## Initialize optional handlers/hit containers to avoid UnboundLocalError when flags not provided
    proteomeEpitopesBestHits = None
    iedbTcruziEpitopesHits = None
    iedbTcruziTableHandler = None
    iedbHumanEpitopesHits = None
    iedbHumanTableHandler = None
    
    #Proteome hits
    if proteome_hits_fp:
        proteomeBlastHandler = BlastResults(proteome_hits_fp)
        logger.info(f"Proteome Parameters: min_length: {ProteomeHistsParams.MIN_LENGTH.value}, min_identity: {ProteomeHistsParams.MIN_IDENTITY.value}, no_gaps: {ProteomeHistsParams.NO_GAPS.value}, qseqid: {ProteomeHistsParams.QSEQID.value}")
        proteomeEpitopesBestHits = proteomeBlastHandler.best_hits_by_qseqid_and_bitscore(ProteomeHistsParams.MIN_LENGTH.value, ProteomeHistsParams.MIN_IDENTITY.value, ProteomeHistsParams.NO_GAPS.value, ProteomeHistsParams.QSEQID.value)
        

    #IEDB options
    validate_iedb_options(iedb, iedb_tcruzi_epitopes_fp, iedb_tcruzi_epitopes_hits_fp, iedb_human_epitopes_hits_fp, iedb_human_epitopes_fp)
    if iedb:
        iedbTcruziBlastHandler = BlastResults(iedb_tcruzi_epitopes_hits_fp)
        logger.info(f"IEDB T. cruzi Parameters: min_length: {DefaultHistsParams.MIN_LENGTH.value}, min_identity: {DefaultHistsParams.MIN_IDENTITY.value}, no_gaps: {DefaultHistsParams.NO_GAPS.value}, qseqid: {DefaultHistsParams.QSEQID.value}")
        iedbTcruziEpitopesHits = iedbTcruziBlastHandler.filter_hits(DefaultHistsParams.MIN_LENGTH.value, DefaultHistsParams.MIN_IDENTITY.value, DefaultHistsParams.NO_GAPS.value, DefaultHistsParams.QSEQID.value)
        iedbTcruziTableHandler = IEDBEpitopeTableHandler(iedb_tcruzi_epitopes_fp)

        ## Human IEDB Epitopes
        if iedb_human_epitopes_hits_fp:
            iedbHumanBlastHandler = BlastResults(iedb_human_epitopes_hits_fp)
            logger.info(f"IEDB Human Parameters: min_length: {DefaultHistsParams.MIN_LENGTH.value}, min_identity: {DefaultHistsParams.MIN_IDENTITY.value}, no_gaps: {DefaultHistsParams.NO_GAPS.value}, qseqid: {DefaultHistsParams.QSEQID.value}")
            iedbHumanEpitopesHits = iedbHumanBlastHandler.filter_hits(DefaultHistsParams.MIN_LENGTH.value, DefaultHistsParams.MIN_IDENTITY.value, DefaultHistsParams.NO_GAPS.value, DefaultHistsParams.QSEQID.value)
            iedbHumanTableHandler = IEDBEpitopeTableHandler(iedb_human_epitopes_fp)
 
    inserts_group_set = load_inserts_group(inserts_group_fp)
    
    
    
    logger.debug(f"CC_idx\tNof_Nodes\tNof_Edges\tNof_Epitopes\tNof_Unique_peptides\tNof_Unique_Reads\tEpitope_candidates\tGenomic_Region_Annotation")

    cc_sorted_list = get_sorted_connected_components(graph)
    epitope_candidates_graph = []
    report_json_epitopes_list = []
    for cc in cc_sorted_list:
        logger.debug(f"Connected component: {cc}")        
        reads_cc = set()
        pepiteds_cc = set()
        genomic_region_annotation_cc = []
        epitope_candidates_cc = []
        genomic_region_locus = []
        for node, attributes in cc.nodes(data=True):
            if not single_reads and len(attributes.get('reads')) < 2:
                # Skipping because not allowed single read on report
                logger.info(f"Skipping because not allowed single read on report: Node {node} - component id: {attributes.get('component_id')}")
                continue

            cc_ids = set()
            genomic_region_locus.append(_get_locus_from_segment(node))
        
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
                    genomic_region_annotation_cc.extend(parse_genomic_region_annotation(value))

            if len(cc_ids) > 1:
                raise Exception(f"CC {cc} has multiple component_ids: {cc_ids}")
            
            cc_id = cc_ids.pop()

            msa_page_name = _create_mview_name(cc_id, pepiteds_cc)

        if single_reads:
            epitope_candidates_graph.extend(epitope_candidates_cc)
        else:
            if len(reads_cc) >= 2: 
                epitope_candidates_graph.extend(epitope_candidates_cc)
            
        
        
        ept_candidantes_seq = []
        for e in epitope_candidates_cc:
            ept_candidantes_seq.append(str(e.seq))
            logger.debug(f"{e.id}\t{cc.number_of_nodes()}\t{cc.number_of_edges()}\t{len(epitope_candidates_cc)}\t{len(pepiteds_cc)}\t{len(reads_cc)}\t{str(e.seq)}\t{genomic_region_annotation_cc}")
            
            
            if single_reads:
                create_cc_json_entity(proteomeEpitopesBestHits, iedb, iedbTcruziEpitopesHits, iedbTcruziTableHandler, iedbHumanEpitopesHits, iedbHumanTableHandler, report_json_epitopes_list, cc, reads_cc, pepiteds_cc, genomic_region_annotation_cc, genomic_region_locus, msa_page_name, e, inserts_group_set)
            else:
                if len(reads_cc) >= 2: 
                    create_cc_json_entity(proteomeEpitopesBestHits, iedb, iedbTcruziEpitopesHits, iedbTcruziTableHandler, iedbHumanEpitopesHits, iedbHumanTableHandler, report_json_epitopes_list, cc, reads_cc, pepiteds_cc, genomic_region_annotation_cc, genomic_region_locus, msa_page_name, e, inserts_group_set)

    
    output_file_name = _get_output_file_name(graph_file, outdir)
    _write_output_fasta_file(epitope_candidates_graph, f"{output_file_name}.fasta")
    
    # Write the JSON object to a file
    with open(f"{output_file_name}.json", "w") as json_file:
        json.dump(report_json_epitopes_list, json_file, indent=4)

def load_inserts_group(inserts_group_fp):
    inserts_group = {}
    if inserts_group_fp:
        with open(inserts_group_fp, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) != 2:
                    logger.warning(f"Invalid line in inserts group file: {line.strip()}")
                    continue
                read_id, group_id = parts
                inserts_group[read_id] = group_id
    return inserts_group

        

def create_cc_json_entity(proteomeEpitopesBestHits, iedb, iedbTcruziEpitopesHits, iedbTcruziTableHandler, iedbHumanEpitopesHits, iedbHumanBlastHandler, report_json_epitopes_list, cc, reads_cc, pepiteds_cc, features_cc, genomic_region_locus, msa_page_name, e, inserts_group_set):
    ## Default values so JSON fields are always present even when optional inputs are missing
    epitope_proteome_best_hit = None
    tcruziEpitopeIEDBHits = []
    humanEpitopeIEDBHits = []

    if proteomeEpitopesBestHits:
        epitope_proteome_best_hit = retrive_epitope_proteome_best_hit(proteomeEpitopesBestHits, e)
        
    
    if iedb:
        tcruziEpitopeIEDBHits = retrieve_iedb_epitope_details(iedbTcruziEpitopesHits, iedbTcruziTableHandler, e)

        if iedbHumanEpitopesHits:
            humanEpitopeIEDBHits = retrieve_iedb_epitope_details(iedbHumanEpitopesHits, iedbHumanBlastHandler, e)    
    
    number_of_inserts_by_group = create_insert_by_group(reads_cc, inserts_group_set)

                
                


                # Create a JSON object with the same information
    json_data = create_json_epitope_data(cc, reads_cc, number_of_inserts_by_group, pepiteds_cc, features_cc, genomic_region_locus, msa_page_name, e, epitope_proteome_best_hit, tcruziEpitopeIEDBHits, humanEpitopeIEDBHits)
    report_json_epitopes_list.append(json_data)

def create_insert_by_group(reads_cc, inserts_group_set):
    group_counts = {}
    if inserts_group_set:
        # Count number of insert_ids by group
        for r in reads_cc:
            insert_id = r.split('_', 1)[0]
            group = inserts_group_set.get(f"@{insert_id}", "NoGroup")
            group_counts[group] = group_counts.get(group, 0) + 1

    return group_counts

def create_json_epitope_data(cc, reads_cc, number_of_inserts_by_group, pepiteds_cc, features_cc, genomic_region_locus, msa_page_name, e, epitope_proteome_best_hit, tcruziEpitopeIEDBHits, humanEpitopeIEDBHits):
    json_data = {
                "ID": e.id,
                "Number of Genomic Regions": len(genomic_region_locus),
                "Number of Peptides": len(pepiteds_cc),  # unique peptides
                "Number of Inserts": len(reads_cc),  # unique reads
                "Number of Inserts by Group": number_of_inserts_by_group,
                "Epitope": str(e.seq),
                "MSA": msa_page_name,
                "Genomic Region Locus": genomic_region_locus,
                "Features": {
                    "ProteinBestHit": epitope_proteome_best_hit,
                    "GenomicRegionsAnnotation": features_cc,
                    "TCruziIEDB": tcruziEpitopeIEDBHits,
                    "HumanIEDB": humanEpitopeIEDBHits,
                }
            }
    
    return json_data

def retrieve_iedb_epitope_details(iedbEpitopesHits, IEDBTableHandler, e):
    iedb_epitopes_info = []
    epitope_pos = BlastColumns.QSEQID.value
    iedb_epitope_pos = BlastColumns.SSEQID.value

    for iedbHits in iedbEpitopesHits:
        if(iedbHits[epitope_pos] == e.id):
            iedbEpitopeInfo = IEDBTableHandler.get_by_epitope_id(iedbHits[iedb_epitope_pos])
            if iedbEpitopeInfo is None:
                logger.warning(f"No IEDB information found for epitope ID: {iedbHits[iedb_epitope_pos]}")
                continue
            logger.debug(f"IEDB Best hit for epitope {iedbHits}")
            logger.debug(f"{iedbEpitopeInfo['Epitope_id']}\t{iedbEpitopeInfo['epitope_sequence']}\t{iedbEpitopeInfo['Epitope - Source Molecule']}\t{iedbEpitopeInfo['Epitope - Source Molecule IRI']}")
            iedb_epitopes_info.append({"IEDB_id": iedbEpitopeInfo['Epitope_id'],
                                              "sequence": iedbEpitopeInfo['epitope_sequence'],
                                              "qstart": iedbHits[BlastColumns.QSTART.value],
                                              "qend": iedbHits[BlastColumns.QEND.value],
                                              "sstart": iedbHits[BlastColumns.SSTART.value],
                                              "send": iedbHits[BlastColumns.SEND.value],
                                              "source_molecule": iedbEpitopeInfo['Epitope - Source Molecule'],
                                              "source_molecule_IRI": iedbEpitopeInfo['Epitope - Source Molecule IRI']})
    return iedb_epitopes_info

def retrive_epitope_proteome_best_hit(proteomeEpitopesBestHits, e):
    epitope_proteome_best_hit_info = {}
    epitope_pos = BlastColumns.QSEQID.value
    protein_pos = BlastColumns.SSEQID.value
    epbh = proteomeEpitopesBestHits.get(e.id)
    if epbh:
        epitope_proteome_best_hit_info = {"epitope_id": epbh[epitope_pos],
                                          "protein_id": epbh[protein_pos],
                                          "protein_description": remove_protein_id_and_taxa_info(epbh[BlastColumns.STITLE.value], epbh[protein_pos]),
                                          "pident": epbh[BlastColumns.PIDENT.value],
                                          "length": epbh[BlastColumns.LENGTH.value],
                                          "mismatches": epbh[BlastColumns.MISMATCH.value],
                                          "gap_opens": epbh[BlastColumns.GAPOPEN.value],
                                          "qstart": epbh[BlastColumns.QSTART.value],
                                          "qend": epbh[BlastColumns.QEND.value],
                                          "sstart": epbh[BlastColumns.SSTART.value],
                                          "send": epbh[BlastColumns.SEND.value],
                                          "evalue": epbh[BlastColumns.EVALUE.value],
                                          "bitscore": epbh[BlastColumns.BITSCORE.value]}
    return epitope_proteome_best_hit_info 

def remove_protein_id_and_taxa_info(originalDescription: str, protein_id: str):
    return originalDescription.replace(protein_id, "").replace("[Trypanosoma cruzi]", "").strip()
    
        
def validate_iedb_options(iedb, iedb_tcruzi_epitopes_fp, iedb_tcruzi_epitopes_hits_fp, iedb_human_epitopes_hits_fp, iedb_human_epitopes_fp):
    if iedb:
        if not iedb_tcruzi_epitopes_fp or not iedb_tcruzi_epitopes_hits_fp:
            logger.error(
                "When --iedb is set, you must also provide both  --iedb-tcruzi-epitopes-hits/-tcruzi-epitopes-hits and --iedb-tcruzi-epitopes/-tcruzi-epitopes  options."                
            )
            raise click.UsageError(
                "Missing required options: --iedb-tcruzi-epitopes/-tcruzi-epitopes and/or --iedb-tcruzi-epitopes-hits/-tcruzi-epitopes-hits when --iedb is used."
            )
  
    if iedb_human_epitopes_hits_fp and not iedb_human_epitopes_fp:
        logger.error(
            "When --iedb-human-epitopes-hits/-human-epitopes-hits is set, you must also provide --iedb-human-epitopes/-human-epitopes option."
        )
        raise click.UsageError(
            "Missing required option: --iedb-human-epitopes/-human-epitopes when --iedb-human-epitopes-hits/-human-epitopes-hits is used."
        )
    
def parse_genomic_region_annotation(genomic_region_annotation_str):
    # Example input: ["WNWZ01000121.1-8261-8334 | KAF8303398.1 | tryptophanyl-tRNA synthetase, putative | 1.0", "WNWZ01000189.1-16368-16457 | pseudogene | tryptophanyl-tRNA synthetase, pseudonote | 0.8640776699029126"]
    genomic_region_annotation = []
    for gra in genomic_region_annotation_str: 
        if gra:
            parts = [p.strip() for p in gra.split('|')]
            if len(parts) != 4:
                raise ValueError(f"Invalid genomic region annotation string format: {gra}")
            genomic_region_annotation.append({
        "genomic_region": parts[0],
        "type": parts[1],
        "description": parts[2],
        "coverage": float(parts[3])
            })
    
    return genomic_region_annotation

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


if __name__ == "__main__":
    main()