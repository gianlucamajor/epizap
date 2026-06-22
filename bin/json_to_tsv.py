import json
import pandas as pd

def converter_json_para_tsv(arquivo_entrada, arquivo_saida):
    # Carregar os dados do arquivo JSON
    with open(arquivo_entrada, 'r', encoding='utf-8') as f:
        dados = json.load(f)

    linhas = []
    
    for item in dados:
        # Extração dos campos principais
        num_peptides = item.get('Number of Peptides', 0)
        
        linha = {
            'ID': item.get('ID', ''),
            'PE': item.get('Epitope', ''),
            'Num_AGR': item.get('Number of Genomic Regions', 0),
            'Num_Peptides': num_peptides,
            'Num_Inserts': item.get('Number of Inserts', 0),
            'MSA': item.get('MSA', '')
        }

        # Classificação automática do tipo de GAGR
        if isinstance(num_peptides, int):
            linha['GAGR_Type'] = 'GAGR-SP' if num_peptides == 1 else ('GAGR-MP' if num_peptides > 1 else '')
        else:
            linha['GAGR_Type'] = ''

        # Achatando o dicionário 'Number of Inserts by Group'
        grupos_insercao = item.get('Number of Inserts by Group', {})
        linha['Inserts_Control'] = grupos_insercao.get('control_a', 0)
        linha['Inserts_Asympto'] = grupos_insercao.get('asympto_a', 0)
        linha['Inserts_CCC_Mild'] = grupos_insercao.get('CCC_mild_a', 0)
        linha['Inserts_CCC_Severe'] = grupos_insercao.get('CCC_severe_a', 0)

        # Achatando os dados do 'ProteinBestHit' dentro de 'Features'
        features = item.get('Features', {})
        best_hit = features.get('ProteinBestHit', {})
        linha['Protein_ID'] = best_hit.get('protein_id', '')
        linha['Protein_Description'] = best_hit.get('protein_description', '')
        linha['BestHit_Evalue'] = best_hit.get('evalue', '')
        linha['BestHit_Pident'] = best_hit.get('pident', '')

        # Convertendo a lista de loci em uma string separada por ponto e vírgula
        locus_lista = item.get('Genomic Region Locus', [])
        linha['Genomic_Region_Locus'] = '; '.join(locus_lista) if locus_lista else ''

        linhas.append(linha)

    # Criar um DataFrame do Pandas e exportar para TSV
    df = pd.DataFrame(linhas)
    df.to_csv(arquivo_saida, sep='\t', index=False, encoding='utf-8')
    print(f"Conversão concluída! Arquivo salvo em: {arquivo_saida}")

# Execução do script
if __name__ == "__main__":
    # Substitua pelos caminhos reais caso os arquivos estejam em outras pastas
    json_input = 'epitopes-data.json'
    tsv_output = 'epitopes-data.tsv'
    
    converter_json_para_tsv("/home/gianluca/Documents/epitopes-data.json", "/home/gianluca/Documents/epitopes-data.tsv")