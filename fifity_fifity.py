#!/usr/bin/env python3
# coding: utf-8
"""
Extrai, para cada gene (Mus musculus), os escores dos 50 primeiros e dos
50 últimos códons do CDS usando uma tabela CSV (códon → valor).

• Respeita rate-limit das E-utilities (delay + retry/back-off).
• Mostra o CDS completo no LOG (console + arquivo rotativo).
• Gera, para cada coluna do Excel:
    - <col>_first50_<timestamp>.csv
    - <col>_last50_<timestamp>.csv
    - <col>_errors_<timestamp>.csv
"""
from __future__ import annotations

import csv
import logging
import logging.handlers
import os
import re
import socket
import time
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Iterable, Tuple

import pandas as pd
from Bio import Entrez, SeqIO
from rich.console import Console
from rich.progress import Progress, TimeElapsedColumn, TimeRemainingColumn
from urllib.error import HTTPError, URLError
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception

# ─────────────── CONFIGURAÇÃO ───────────────
BASE       = Path(__file__).parent
DATA_FILE  = BASE / "data" / "upregulated_proteins.xlsx"
TABLE_FILE = BASE / "data" / "tabela.csv"

OUT_DIR = BASE / "product_codons";  OUT_DIR.mkdir(exist_ok=True)
LOG_DIR = BASE / "LOG_codons";      LOG_DIR.mkdir(exist_ok=True)

STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

Entrez.email   = "seu_email@exemplo.com"  # substitua
Entrez.api_key = ""                       # opcional
NCBI_DELAY     = 0.12 if Entrez.api_key else 0.34
socket.setdefaulttimeout(30)

# ─────────────── LOGGING ─────────────────────
console = Console()
# Console colorido + arquivo rotativo
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%H:%M:%S"))
fh = logging.handlers.RotatingFileHandler(LOG_DIR / "run.log", maxBytes=5_000_000, backupCount=3, encoding="utf-8")
fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-7s | %(name)s | %(message)s", "%Y-%m-%d %H:%M:%S"))

logging.basicConfig(level=logging.INFO, handlers=[sh, fh])
log = logging.getLogger("codon50")

@contextmanager
def timed(msg: str):
    log.info(f"{msg}…")
    t0 = time.perf_counter()
    yield
    log.info(f"{msg} concluído em {time.perf_counter() - t0:.1f}s")

# ─────────────── TABELA DE CÓDONS ─────────────────────
def load_codon_table(path: Path) -> dict[str, float]:
    df = pd.read_csv(path, sep=";", decimal=",", engine="python")
    cod_col, val_col = df.columns[:2]
    tbl = {row[cod_col].strip().upper(): row[val_col] for _, row in df.iterrows()}
    log.info(f"Tabela {path.name} carregada ({len(tbl)} códons)")
    return tbl

CODON_SCORE = load_codon_table(TABLE_FILE)

# ─────────────── UTILITÁRIOS ─────────────────────
def snake_case(s: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^0-9a-zA-Z]+", "_", s)).strip("_").lower()

def polite_call(fn, *args, **kw):
    h = fn(*args, **kw)
    time.sleep(NCBI_DELAY)
    return h

def is_retryable(exc: Exception) -> bool:
    from requests.exceptions import RequestException
    return isinstance(exc, (HTTPError, URLError, RequestException, TimeoutError, socket.timeout))

# ─────────────── EXCEL → LISTAS ─────────────────────
def load_lists() -> dict[str, list[str]]:
    df = pd.read_excel(DATA_FILE, sheet_name=0)
    df.columns = map(snake_case, df.columns)
    for col in df.select_dtypes("object"):
        df[col] = (
            df[col]
            .fillna("")
            .astype(str)
            .str.split(r"[;,]")
            .apply(lambda lst: [g.strip() for g in lst if g.strip()])
        )
    return {c: [g for sub in df[c] for g in sub] for c in df.select_dtypes("object")}

# ─────────────── NCBI HELPERS ─────────────────────
TERM_TEMPLATES = (
    "Mus musculus[ORGN] AND {gene}[Gene] AND transcript variant 1[Title] AND mRNA[Filter]",
    "Mus musculus[ORGN] AND {gene}[Gene] AND mRNA[Filter]",
)

@retry(stop=stop_after_attempt(5),
       wait=wait_exponential(multiplier=1, min=2, max=32),
       retry=retry_if_exception(is_retryable))
def find_uid(gene: str) -> str | None:
    for tmpl in TERM_TEMPLATES:
        handle = polite_call(Entrez.esearch,
                             db="nucleotide",
                             term=tmpl.format(gene=gene),
                             retmode="xml",
                             retmax=1)
        rec = Entrez.read(handle); handle.close()
        ids = rec.get("IdList", [])
        if ids:
            return ids[0]
    return None

@retry(stop=stop_after_attempt(5),
       wait=wait_exponential(multiplier=1, min=2, max=32),
       retry=retry_if_exception(is_retryable))
def fetch_gb(uid: str) -> SeqIO.SeqRecord:
    h = polite_call(Entrez.efetch,
                    db="nucleotide", id=uid,
                    rettype="gb", retmode="text")
    rec = SeqIO.read(h, "genbank"); h.close()
    return rec

def cds_bounds(rec: SeqIO.SeqRecord) -> Tuple[int, int] | Tuple[None, None]:
    for feat in rec.features:
        if feat.type == "CDS":
            return int(feat.location.start), int(feat.location.end)
    return None, None

# ─────────────── PROCESSAMENTO DE CADA GENE ─────────────────────
def process_gene(gene: str) -> Tuple[dict, str] | Tuple[None, str]:
    # 1) UID
    uid = find_uid(gene)
    if uid is None:
        return None, "UID não encontrado"
    # 2) Fetch GenBank
    rec = fetch_gb(uid)
    seq = str(rec.seq).upper()
    # 3) Bounds do CDS
    s, e = cds_bounds(rec)
    if s is None:
        return None, "CDS não anotado"
    # 4) Extrai CDS e log completo
    cds_seq = seq[s:e]
    log.info(f"[{gene}] CDS completo ({len(cds_seq)} nt): {cds_seq}")
    # 5) Lista de escores
    codons = [cds_seq[i:i+3] for i in range(0, len(cds_seq)-2, 3)]
    vals   = [CODON_SCORE[c] for c in codons if c in CODON_SCORE]
    if not vals:
        return None, "Nenhum códon pontuado"
    # 6) Separa primeiros 50 / últimos 50
    first_vals = vals[:50]
    last_vals  = vals[-50:] if len(vals) >= 50 else vals
    # Empacota para string CSV
    pack = lambda lst: ";".join(f"{v:.6g}" for v in lst)
    return {
        "gene": gene,
        "first50": pack(first_vals),
        "last50":  pack(last_vals),
    }, None

# ─────────────── AJUDA CSV ─────────────────────
FIELDS_FIRST = ["gene", "first50"]
FIELDS_LAST  = ["gene", "last50"]
def write_csv(path: Path, rows: list[dict], fields: list[str]):
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    log.info(f"Wrote {path.name} ({len(rows)} registros)")

# ─────────────── MAIN ─────────────────────
def main():
    with timed("Carregando listas do Excel"):
        lists = load_lists()

    for col, genes in lists.items():
        base = snake_case(col)
        log.info(f"→ Processando coluna '{col}' ({len(genes)} genes)")
        first_rows, last_rows, errs = [], [], []

        with Progress(
            "[progress.percentage]{task.percentage:>3.0f}%",
            "•", "{task.description}",
            TimeElapsedColumn(), "/", TimeRemainingColumn(),
            console=console, transient=True
        ) as prog:
            task = prog.add_task(f"[cyan]{base}", total=len(genes))
            for gene in genes:
                try:
                    res, err = process_gene(gene)
                except Exception as ex:
                    res, err = None, f"erro inesperado: {ex}"
                if res:
                    first_rows.append({"gene": res["gene"], "first50": res["first50"]})
                    last_rows.append( {"gene": res["gene"], "last50":  res["last50"]})
                else:
                    errs.append({"gene": gene, "error": err})
                prog.advance(task)

        # Grava dois CSVs + erros
        write_csv(OUT_DIR / f"{base}_first50_{STAMP}.csv", first_rows, FIELDS_FIRST)
        write_csv(OUT_DIR / f"{base}_last50_{STAMP}.csv",  last_rows,  FIELDS_LAST)
        write_csv(LOG_DIR / f"{base}_errors_{STAMP}.csv", errs, ["gene", "error"])

    log.info("🏁 Todos os genes foram processados com sucesso.")

if __name__ == "__main__":
    main()
