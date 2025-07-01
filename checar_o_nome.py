#!/usr/bin/env python3
# coding: utf-8
"""
Calcula, para cada gene (Mus musculus), a posição do CDS, a sequência Kozak
e a média geométrica dos códons (tabela.csv).
Prefere “transcript variant 1”; se não encontrar, cai para qualquer mRNA.

• Respeita o rate‑limit das E‑utilities (delay automático + back‑off).
• Usa sua API‑key (até 10 req/s).
• Gera um CSV por coluna do Excel + um CSV de erros.
• Logging detalhado no console (cores + barra de progresso) e em arquivo rotativo.
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
from statistics import geometric_mean
from typing import Iterable, Tuple

import pandas as pd
from Bio import Entrez, SeqIO
from rich.console import Console
from rich.progress import Progress, TimeElapsedColumn, TimeRemainingColumn
from urllib.error import HTTPError, URLError
from tenacity import (
    retry,
    stop_after_attempt,
    wait_exponential,
    retry_if_exception,
)

# ─────────────────────────── Configuração global ────────────────────────────
BASE        = Path(__file__).parent
DATA_FILE   = BASE / "data" / "upregulated_proteins.xlsx"
TABLE_FILE  = BASE / "data" / "tabela.csv"

OUT_DIR     = BASE / "product_kozak"
ERR_DIR     = BASE / "LOG_kozak"
LOG_DIR     = ERR_DIR                         # mesmo diretório dos erros

OUT_DIR.mkdir(exist_ok=True)
ERR_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)

STAMP       = datetime.now().strftime("%Y%m%d_%H%M%S")

Entrez.email   = "seu_email@exemplo.com"                      # ajuste!
Entrez.api_key = "f28d31a184e32070782a50e14e525008de07"       # sua chave

# delay mínimo entre chamadas (API‑key ⇒ 10 req/s; sem chave ⇒ 3 req/s)
NCBI_DELAY = 0.12 if Entrez.api_key else 0.34

socket.setdefaulttimeout(30)                                  # timeout de leitura

# ─────────────────────────── Logging ────────────────────────────
def setup_logging(level: str = "INFO") -> Console:
    console = Console()
    level   = getattr(logging, level.upper(), logging.INFO)

    # Handler colorido no console
    sh = logging.StreamHandler()
    sh.setLevel(level)
    sh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s",
                                      "%H:%M:%S"))

    # Handler de arquivo rotativo
    fh = logging.handlers.RotatingFileHandler(
        LOG_DIR / "run.log", maxBytes=2_000_000, backupCount=3, encoding="utf-8"
    )
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
        "%Y-%m-%d %H:%M:%S"
    ))

    logging.basicConfig(level=level, handlers=[sh, fh], force=True)
    console.log(f"[bold green]Log level[/] → [yellow]{logging.getLevelName(level)}[/]")
    return console


console = setup_logging(os.getenv("LOG_LEVEL", "INFO"))
log     = logging.getLogger("kozak")

@contextmanager
def timed(msg: str):
    log.info(msg + "…")
    t0 = time.perf_counter()
    try:
        yield
    finally:
        log.info("%s concluído em %.1f s", msg, time.perf_counter() - t0)

# ─────────────────────────── Tabela de códons ───────────────────────────────
def load_codon_table(path: Path) -> dict[str, float]:
    df = pd.read_csv(path, sep=";", decimal=",", engine="python")
    codon_col, val_col = df.columns[:2]
    tbl = {row[codon_col].strip().upper(): row[val_col] for _, row in df.iterrows()}
    log.info("Loaded %d codon scores from %s", len(tbl), path.name)
    return tbl


CODON_SCORE = load_codon_table(TABLE_FILE)

# ─────────────────────────── Utilitários ─────────────────────────────────────
def snake_case(s: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^0-9a-zA-Z]+", "_", s)).strip("_").lower()


def polite_call(fn, *args, **kw):
    """Executa E‑utility com pausa respeitando o limite de taxa."""
    handle = fn(*args, **kw)
    time.sleep(NCBI_DELAY)
    return handle


def geom_mean(vals: Iterable[float]) -> float | None:
    lst = list(vals)
    return geometric_mean(lst) if lst else None


def is_retryable(exc: Exception) -> bool:
    from requests.exceptions import RequestException
    return isinstance(
        exc,
        (HTTPError, URLError, RequestException, TimeoutError, socket.timeout),
    )

# ─────────────────────────── Excel → listas de genes ────────────────────────
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

# ─────────────────────────── NCBI helpers ───────────────────────────────────
TERM_TEMPLATES = (
    "Mus musculus[ORGN] AND {gene}[Gene] AND transcript variant 1[Title] AND mRNA[Filter]",
    "Mus musculus[ORGN] AND {gene}[Gene] AND mRNA[Filter]",
)

@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=2, max=32),
    retry=retry_if_exception(is_retryable),
)
def find_uid(gene: str) -> Tuple[str | None, str]:
    """Retorna (UID, estratégia) ou (None, motivo)."""
    for i, tmpl in enumerate(TERM_TEMPLATES, 1):
        term = tmpl.format(gene=gene)
        h    = polite_call(Entrez.esearch,
                           db="nucleotide", term=term, retmode="xml", retmax=1)
        rec  = Entrez.read(h); h.close()
        ids  = rec.get("IdList", [])
        if ids:
            return ids[0], f"strategy_{i}"
    return None, "no_uid"


@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=2, max=32),
    retry=retry_if_exception(is_retryable),
)
def fetch_gb(uid: str) -> SeqIO.SeqRecord:
    h = polite_call(Entrez.efetch,
                    db="nucleotide", id=uid, rettype="gb", retmode="text")
    rec = SeqIO.read(h, "genbank"); h.close()
    return rec


def cds_bounds(rec: SeqIO.SeqRecord) -> tuple[int | None, int | None]:
    for feat in rec.features:
        if feat.type == "CDS":
            return int(feat.location.start), int(feat.location.end)
    return None, None

# ─────────────────────────── Processamento de gene ──────────────────────────
def process_gene(gene: str) -> Tuple[dict | None, str | None]:
    uid, strat = find_uid(gene)
    if uid is None:
        return None, f"UID not found ({strat})"

    rec     = fetch_gb(uid)
    seq     = str(rec.seq).upper()
    cds_s, cds_e = cds_bounds(rec)
    if cds_s is None:
        return None, "CDS not annotated"

    cds_seq   = seq[cds_s:cds_e]
    kozak_seq = seq[max(0, cds_s - 5): cds_s + 4]           # –5 … +3 (9 nt)

    codon_vals = [
        CODON_SCORE[codon]
        for codon in (cds_seq[i:i + 3] for i in range(0, len(cds_seq) - 2, 3))
        if codon in CODON_SCORE
    ]
    gm = geom_mean(codon_vals)
    if gm is None:
        return None, "No scored codons"

    return {
        "gene": gene,
        "cds_start": cds_s,
        "geom_mean_cds": round(gm, 4),
        "kozak_seq": kozak_seq,
    }, None

# ─────────────────────────── CSV helpers ────────────────────────────────────
CSV_COLS = ["gene", "cds_start", "geom_mean_cds", "kozak_seq"]


def write_csv(path: Path, rows: list[dict]):
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLS)
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote %s (%d records)", path.name, len(rows))


def write_errors(path: Path, errs: list[dict]):
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "error"])
        w.writeheader()
        w.writerows(errs)
    if errs:
        log.warning("Wrote %s (%d errors)", path.name, len(errs))

# ─────────────────────────── Main ───────────────────────────────────────────
def main():
    t0 = time.perf_counter()
    with timed("Carregando listas do Excel"):
        lists = load_lists()

    for col, genes in lists.items():
        ok_rows: list[dict] = []
        err_rows: list[dict] = []
        base = snake_case(col)

        log.info("→ Processando lista '%s' (%d genes)", col, len(genes))
        with Progress(
            "[progress.percentage]{task.percentage:>3.0f}%",
            "•",
            "{task.description}",
            TimeElapsedColumn(),
            "/",
            TimeRemainingColumn(),
            console=console,
            transient=True,
        ) as progress:
            task = progress.add_task(f"[cyan]{base}", total=len(genes))
            for g in genes:
                try:
                    res, err = process_gene(g)
                except Exception as exc:            # captura falha inesperada
                    err = f"unexpected: {exc}"
                    res = None
                if res:
                    ok_rows.append(res)
                    log.debug("%s → OK (gm=%s)", g, res["geom_mean_cds"])
                else:
                    err_rows.append({"gene": g, "error": err})
                    log.warning("%s → %s", g, err)
                progress.advance(task)

        write_csv (OUT_DIR / f"{base}_{STAMP}.csv", ok_rows)
        write_errors(ERR_DIR / f"{base}_errors_{STAMP}.csv", err_rows)
        log.info("✔ %s: %d ok, %d erros\n", col, len(ok_rows), len(err_rows))

    log.info("🏁 Completo em %.1fs", time.perf_counter() - t0)


if __name__ == "__main__":
    main()
