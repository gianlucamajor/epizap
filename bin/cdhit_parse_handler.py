import click
import logging.handlers

from core.CDHitParser import CDHitParser

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger("CDHIT-PARSER")
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

@click.command(help="Aim of this program is decrease redundancy of the peptides (epitopes) in the same Graph Connected Component (CC).")
@click.argument("files", nargs=-1, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--peptides", "-pep", type=click.Path(), required=True, help="The path to the fasta file with the peptides.")
@click.option("--outdir", "-o", type=click.Path(), default="results_graph_21_01_2025/epitopes-by-msa-core-and-lonely/")
def main(files:click.Path, peptides:click.Path, outdir:click.Path):

    epitopes_sequences_no_clustered = []
    epitopes_sequences = read_fasta(peptides)
    
    
    for f in files:
        parser = CDHitParser(f)
        clusters =  parser.get_clusters_with_more_than_one_peptides()
        peptides_id_no_cluster = parser.get_peptides_no_cluster()    
        logger.info(f"CC ID: {parser.cc_id}")
        logger.info(f"NUMBER IDS NO CLUSTER: {parser.get_total_peptides_no_cluster()}, NUMBER OF CLUSTERS {len(clusters)}")
        logger.debug(f"IDS NO CLUSTER: {peptides_id_no_cluster}")
        logger.debug(f"CLUSTERS: {clusters}")
        [logger.debug(f"Cluster {cluster_id}: number of peptides:{len(sequences)}") for cluster_id, sequences in clusters.items()]
        

        clusters_id = list(clusters.keys())
        for cluster_id in clusters_id:
            epitopes_sequences_clustered = []
            for peptide_id in clusters[cluster_id]:
                epitopes_sequences_clustered.append(epitopes_sequences[peptide_id])
            if len(epitopes_sequences_clustered) > 0:
                write_output_file_(epitopes_sequences_clustered, f"{outdir}/cc-clustered/epitopes-from-cc-{parser.cc_id}-cluster-{cluster_id}.fasta")
      

        for peptide_id in peptides_id_no_cluster:
            epitopes_sequences_no_clustered.append(epitopes_sequences[peptide_id])

        logger.info("_"*100)
    
    write_output_file_(epitopes_sequences_no_clustered, outdir+"/epitopes-on-same-cc-no-cluster.fasta")

    
    


def write_output_file_(sequences, output_file_name):
    with open(output_file_name, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

def read_fasta(faa_file):
    return SeqIO.index(faa_file, "fasta")

if __name__ == "__main__":
    main()
        

