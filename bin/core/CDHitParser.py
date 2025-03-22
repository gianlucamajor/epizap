#!/usr/bin/env python
import re

class CDHitParser:
    def __init__(self, filename):
        self.filename = filename
        self.cdhit_clusters = {}
        self.cc_id = None
        
        self.parse_cdhit_clusters()
        self._get_cc_id()


    def parse_cdhit_clusters(self):
        current_cluster = None
        with open(self.filename, 'r') as f:
            for line in f:
                if line.startswith(">Cluster"):
                    current_cluster = int(line.split()[1])  # Extract cluster number
                    self.cdhit_clusters[current_cluster] = []
                else:
                    if current_cluster is not None:
                        self.cdhit_clusters[current_cluster].append(line.strip())  # Add sequence ID
    

    def get_cdhit_clusters(self):
        return self.cdhit_clusters
    
    def get_number_of_cdhit_clusters(self):
        return len(self.clusters)
    
    def _extract_seq_id(self, sequence):
        # Use regex to find the sequence ID
        # line examples:
        # "9       15aa, >7352-20... at 100.00%",
        # "10      49aa, >major_1... *"
        match = re.search(r'>([\w-]+)', sequence)
        if match:
            return match.group(1)  # Return the captured sequence ID
        return None
    
    # return all alone peptides in own cluster CD-HIT cluster
    def get_peptides_no_cluster(self):
        peptides_id_no_cluster = []
        for cluster_id, sequences in self.cdhit_clusters.items():
            if len(sequences) == 1: 
                peptides_id = self._extract_seq_id(sequences[0])
                peptides_id_no_cluster.append(peptides_id)
        return peptides_id_no_cluster
    
    def get_clusters_with_more_than_one_peptides(self):
        peptides_id_clustering = {}
        for cluster_id, sequences in self.cdhit_clusters.items():
            if len(sequences) > 1:
                peptide_ids_on_cluster = []
                for sequence in sequences:
                    peptide_ids_on_cluster.append(self._extract_seq_id(sequence))
                peptides_id_clustering[cluster_id] = peptide_ids_on_cluster

        return peptides_id_clustering
    
    def get_total_peptides_no_cluster(self):
        return len(self.get_peptides_no_cluster())
    
    def _get_cc_id(self):
        # Use regex to find the Conex Componet ID
        match = re.search(r'epitopes-cc-(\d+)', self.filename)
        if match:
            self.cc_id = match.group(1)  # Return the captured cluster ID

    def print_clusters(self):
        for cluster_id, sequences in self.cdhit_clusters.items():
            print(f"Cluster {cluster_id}:")
            for sequence in sequences:
                print(f"  {sequence}")