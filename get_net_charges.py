#!/usr/bin/env python3
# coding: utf-8
"""
Script integrado com timeouts e delays:
1. Carrega e limpa o Excel (listas por coluna).
2. Busca sequências no NCBI com timeout de socket e delay, evita travamentos.
3. Calcula cargas e janelas deslizantes.
4. Gera CSVs de 30 resíduos, 100 janelas e erros por coluna.
"""
import re
import csv
import logging
import socket
import time
from pathlib import Path
from datetime import datetime

import pandas as pd
from Bio import Entrez, SeqIO
from tenacity import retry, stop_after_attempt, wait_fixed

# ------------------------------------
# Configuração geral
# ------------------------------------
# Timeout global para chamadas de rede (em segundos)
socket.setdefaulttimeout(10)

Entrez.email = "seu_email@exemplo.com"  # substitua pelo seu e-mail
BASE_DIR = Path(__file__).parent
DATA_FILE = BASE_DIR / 'data' / 'upregulated_proteins.xlsx'
PRODUCT_DIR = BASE_DIR / 'product'
LOG_DIR     = BASE_DIR / 'LOG'
TIMESTAMP   = datetime.now().strftime('%Y%m%d_%H%M%S')

PRODUCT_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)

# Logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Delay entre chamadas ao NCBI (segundos)
NCBI_DELAY = 0.4

# ------------------------------------
# Funções utilitárias
# ------------------------------------

def snake_case(s: str) -> str:
    return re.sub(r'_+', '_', re.sub(r'[^0-9a-zA-Z]+', '_', s)).strip('_').lower()


def load_and_clean(fp: Path) -> pd.DataFrame:
    df = pd.read_excel(fp, sheet_name=0)
    df.columns = [snake_case(c) for c in df.columns]
    for col in df.select_dtypes(include='object'):
        df[col] = (
            df[col]
              .fillna("")
              .astype(str)
              .str.split(r'[;,]')
              .apply(lambda lst: [g.strip() for g in lst if g.strip()])
        )
    return df


def flatten_columns(df: pd.DataFrame) -> dict:
    return {col: [gene for sub in df[col] for gene in sub]
            for col in df.select_dtypes(include='object')}

# ------------------------------------
# Consultas NCBI usando Entrez (XML)
# ------------------------------------
@retry(stop=stop_after_attempt(3), wait=wait_fixed(2))
def search_protein_for_gene(gene: str, organism="mus musculus") -> str:
    term = f"{organism}[ORGN] AND {gene}[Gene]"
    logger.info(f"Searching UID for {gene}")
    handle = Entrez.esearch(db="protein", term=term, retmax=1, retmode="xml")
    rec = Entrez.read(handle)
    handle.close()
    ids = rec.get("IdList", [])
    time.sleep(NCBI_DELAY)
    if not ids:
        logger.warning(f"No UID found for {gene}")
        return None
    return ids[0]

@retry(stop=stop_after_attempt(3), wait=wait_fixed(2))
def get_accession_from_uid(uid: str) -> str:
    logger.info(f"Getting accession for UID {uid}")
    handle = Entrez.esummary(db="protein", id=uid, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    time.sleep(NCBI_DELAY)
    if isinstance(records, list) and records:
        rec0 = records[0]
        return rec0.get('AccessionVersion') or rec0.get('Caption')
    logger.warning(f"No summary record for UID {uid}")
    return None

@retry(stop=stop_after_attempt(3), wait=wait_fixed(2))
def fetch_sequence(accession: str) -> tuple[str, str]:
    logger.info(f"Fetching sequence for {accession}")
    handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
    rec = SeqIO.read(handle, "fasta")
    handle.close()
    time.sleep(NCBI_DELAY)
    return str(rec.seq), rec.description

# ------------------------------------
# Cálculo de cargas e janelas
# ------------------------------------

def compute_charges(seq: str) -> list[int]:
    cmap = {'K':1, 'R':1, 'D':-1, 'E':-1}
    charges = [cmap.get(aa.upper(), 0) for aa in seq]
    if charges: charges[0] = 1
    if len(charges) > 1: charges[-1] = -1
    return charges


def sliding_window_sums(charges: list[int], size: int = 30) -> list[int]:
    return [sum(charges[i:i+size]) for i in range(len(charges) - size + 1)]

# ------------------------------------
# Classe de processamento de proteína
# ------------------------------------
class ProteinProcessor:
    def __init__(self, gene: str):
        self.gene = gene
        self.seq_length = 0
        self.charges = []
        self.windows = []

    def process(self) -> bool:
        uid = search_protein_for_gene(self.gene)
        if not uid:
            return False
        acc = get_accession_from_uid(uid)
        if not acc:
            return False
        seq, _ = fetch_sequence(acc)
        self.seq_length = len(seq)
        self.charges = compute_charges(seq)
        self.windows = sliding_window_sums(self.charges)
        return True

    def row_30(self) -> list:
        return [self.gene, self.seq_length] + self.charges[:30] + [""]*(30 - len(self.charges))

    def row_100(self) -> list:
        return [self.gene, self.seq_length] + self.windows[:100] + [""]*(100 - len(self.windows))

# ------------------------------------
# Escrita de CSV
# ------------------------------------

def write_csv(path: Path, header: list, rows: list) -> None:
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)
    logger.info(f"CSV written: {path}")

# ------------------------------------
# Main
# ------------------------------------
def main():
    df = load_and_clean(DATA_FILE)
    lists = flatten_columns(df)

    for col, genes in lists.items():
        base = snake_case(col)
        rows30, rows100, errs = [], [], []
        for gene in genes:
            proc = ProteinProcessor(gene)
            if not proc.process():
                errs.append([gene, 'not found'])
            else:
                rows30.append(proc.row_30())
                rows100.append(proc.row_100())
        header30 = ['gene', 'seq_length'] + [f'R{i+1}' for i in range(30)]
        header100 = ['gene', 'seq_length'] + [f'W{i+1}' for i in range(100)]
        p30 = PRODUCT_DIR / f"{base}_30residues_{TIMESTAMP}.csv"
        p100 = PRODUCT_DIR / f"{base}_100windows_{TIMESTAMP}.csv"
        perr = LOG_DIR / f"{base}_errors_{TIMESTAMP}.csv"
        write_csv(p30, header30, rows30)
        write_csv(p100, header100, rows100)
        write_csv(perr, ['gene', 'error'], errs)

if __name__ == '__main__':
    main()