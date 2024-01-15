#!/usr/bin/env python
import os


import click
from multiprocessing.pool import Pool
from pathlib import Path
from Bio import SeqIO

@click.command(help="Aim of this program is generate a fasta file with lonely peptide sequences.")
@click.argument("files", nargs=-1, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--outdir", "-o", type=click.Path(), default="lonely/lonely_peptides.fasta")
@click.option("--threads", "-t", type=int, default=1)
def main(files:click.Path, outdir, threads:int):
        
    sequences_renamed = []
    with Pool(threads) as pool:
        for record_renamed in pool.map(_rename_record_id, files):
            sequences_renamed.append(record_renamed)
    _write_output_file(sequences_renamed, outdir)
        
        
def _rename_record_id(file):
    file_name = os.path.splitext(os.path.basename(file))[0]
    record = _read_lonely_fasta(file)
    record.id = "_".join([file_name, record.id])
    return record

def _write_output_file(sequences, output_file_name):
    with open(output_file_name, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

def _read_lonely_fasta(faa_file):
    return SeqIO.read(faa_file, "fasta")

if __name__ == "__main__":
    main()