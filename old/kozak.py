import os
import csv
import logging
from typing import Optional, Tuple
import pandas as pd
from Bio import Entrez, SeqIO

# Configuração do e-mail para acesso ao NCBI
Entrez.email = "seu_email@exemplo.com"  # Substitua pelo seu e-mail

# Diretórios de saída
OUTPUT_DIR = "../saida_kozak"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Configuração do log
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(os.path.join(OUTPUT_DIR, "kozak_analysis.log")),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("kozak")

def fetch_mrna_sequence(gene_name: str, organism: str = "Mus musculus") -> Optional[Tuple[str, str]]:
    """
    Busca a sequência de mRNA para um gene específico no NCBI.
    Retorna a sequência e a descrição, ou None se não for encontrada.
    """
    try:
        search_term = f"{gene_name}[Gene Name] AND {organism}[Organism] AND mRNA[Filter]"
        handle = Entrez.esearch(db="nucleotide", term=search_term, retmode="xml", retmax=1)
        record = Entrez.read(handle)
        handle.close()
        id_list = record.get("IdList", [])
        if not id_list:
            logger.warning(f"Nenhuma sequência encontrada para o gene '{gene_name}'.")
            return None
        seq_id = id_list[0]
        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
        seq_record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(seq_record.seq), seq_record.description
    except Exception as e:
        logger.error(f"Erro ao buscar a sequência para o gene '{gene_name}': {e}")
        return None

def find_kozak_sequence(mrna_seq: str) -> Optional[Tuple[int, str]]:
    """
    Procura pela sequência de Kozak (AUGG) na sequência de mRNA.
    Retorna a posição inicial e a sequência anterior, se encontrada.
    """
    mrna_seq = mrna_seq.upper()
    for i in range(len(mrna_seq) - 3):
        codon = mrna_seq[i:i+4]
        if codon.startswith("AUG") and codon[3] == "G":
            upstream_seq = mrna_seq[max(0, i-10):i]
            return i, upstream_seq
    return None

def process_genes_from_excel(excel_path: str, organism: str = "Mus musculus") -> None:
    """
    Processa uma planilha Excel contendo nomes de genes para identificar a sequência de Kozak.
    Gera um arquivo CSV com os resultados.
    """
    try:
        df = pd.read_excel(excel_path)
    except Exception as e:
        logger.error(f"Erro ao ler o arquivo Excel: {e}")
        return

    if df.empty or df.columns.empty:
        logger.error("A planilha está vazia ou não possui colunas.")
        return

    gene_column = df.columns[0]
    gene_list = df[gene_column].dropna().astype(str).tolist()

    output_file = os.path.join(OUTPUT_DIR, "kozak_results.csv")
    with open(output_file, mode="w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Gene", "Descrição", "Posição Kozak", "Sequência Anterior"])
        for gene in gene_list:
            logger.info(f"Processando gene: {gene}")
            result = fetch_mrna_sequence(gene, organism)
            if result is None:
                writer.writerow([gene, "Sequência não encontrada", "", ""])
                continue
            seq, desc = result
            logger.info(f"Sequência de mRNA para {gene}: {seq}")
            kozak_result = find_kozak_sequence(seq)
            if kozak_result is None:
                writer.writerow([gene, desc, "Kozak não encontrada", ""])
            else:
                pos, upstream = kozak_result
                writer.writerow([gene, desc, pos, upstream])
    logger.info(f"Processamento concluído. Resultados salvos em '{output_file}'.")

if __name__ == "__main__":
    # Substitua pelo caminho do seu arquivo Excel
    excel_file_path = "../data/Upregulated_proteins_BKP.xlsx"
    process_genes_from_excel(excel_file_path)
