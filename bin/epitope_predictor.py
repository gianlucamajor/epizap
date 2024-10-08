#!/usr/bin/env python

import click
import re

@click.command(help="Aim of this program is predict epitope from the multiple sequence align (MSA) core. ")
@click.argument("files", nargs=-1, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--outdir", "-o", type=click.Path(), default="")
def main(files:click.Path, outdir:click.Path):
    consensus_pattern = re.compile(r'^consensus/(\d{2,3})%')
    epitope_pattern = re.compile(r'[A-Z]{5,}')
    
    for f in files:
        # print(f)
        epitopes = get_concensus(f, consensus_pattern, epitope_pattern)
        for k in epitopes:
            print(f"{k}")
            for e in epitopes[k]:
                print(e)

def get_concensus(consensus_file, consensus_pattern, epitope_pattern):
    with open(consensus_file, 'r') as cfile:
        epitopes = {}
        consensus_list = []
        for line in cfile:
            ls = line.strip()
            consensus_found = consensus_pattern.search(ls)
            if consensus_found:
                # print(f"Linha correspondente: {ls}")
                # print(f"Valor encontrado: {consensus_found.group(1)}")
                consensus_value = consensus_found.group(1)
                
                epitopes_found = epitope_pattern.findall(ls)
                if epitopes_found:
                    consensus_list.append(f"{consensus_value}:{epitopes_found}")
        
        epitopes[cfile.name] = consensus_list
        return epitopes    


                

                


if __name__ == "__main__":
    main()