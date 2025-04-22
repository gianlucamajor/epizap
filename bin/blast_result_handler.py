
import click
from Bio import SeqIO

from model.BlastResult import BlastResult

@click.command(help="Aim of this program is ")
@click.argument("input", type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.argument("peptide_file", type=click.Path(exists=True, file_okay=True, dir_okay=False))
def main(input:click.Path, peptide_file:click.Path):
    peptides = _read_fasta(peptide_file)
    br_hits =  load_blast_results(input, peptides)
    print(len(br_hits))

    all_hits_ids = set()
    br_hits = {k: v for k, v in sorted(br_hits.items(), key=lambda item: len(item[1]), reverse=True)}
    for k, v in br_hits.items():
        print(k, peptides.get(k).seq)
        print(len(v))
        for br in v:
            all_hits_ids.add(br.sseqid)
            print(br)
            print("-" * 50)
    print("All hits ids", len(all_hits_ids))

    # print(len(br_hits))
    # s = peptides.get("173-1").seq
    # print(len(s), s)
    # l = br_hits.get("173-1")
    # if l is None:
    #     print("No hits found")
    # else:
    #     print(len(l))
    #     for br in l:
    #         print(br)
    #         print("-" * 50)
    

def load_blast_results(tsv_file, peptides):
    blast_results = {}
    with open(tsv_file, "r") as file:
        for line in file:
            # Create a BlastResult object for each line
            blast_hit = BlastResult(line)
            original_sequence = peptides.get(blast_hit.qseqid)
            if original_sequence is None:
                ValueError(f"Peptide {blast_hit.qseqid} not found in the peptide file.")
            
            # it's not enough to filter. Also is necessary to check length of the original sequence;
            if blast_hit.qseqid != blast_hit.sseqid and blast_hit.pident == 100 and  blast_hit.qcovs == 100 and blast_hit.length == len(original_sequence.seq): 
                qseq_hits = blast_results.get(blast_hit.qseqid, None)
                if qseq_hits is None:
                    new_qseq_hit_list = []
                    new_qseq_hit_list.append(blast_hit)
                    blast_results[blast_hit.qseqid] = new_qseq_hit_list
                else:
                    qseq_hits.append(blast_hit)
                    blast_results[blast_hit.qseqid] = qseq_hits
            
    return blast_results    

def _read_fasta(faa_file):
    return SeqIO.index(faa_file, "fasta")
if __name__ == "__main__":
    main()
