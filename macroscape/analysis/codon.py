import logging
import re
from statistics import geometric_mean
from pathlib import Path

import pandas as pd
from Bio.SeqRecord import SeqRecord

from .. import settings

log = logging.getLogger("codon")


class KozakAnalyzer:
    """Calcula posição do CDS, sequência Kozak e média geométrica de códons."""

    def __init__(self, codon_table: Path | None = None):
        path = Path(codon_table or settings.CODON_TABLE)
        df = pd.read_csv(path, sep=";", decimal=",", engine="python")
        codon_col, val_col = df.columns[:2]
        self._score = {r[codon_col].strip().upper(): r[val_col] for _, r in df.iterrows()}
        log.info("Tabela de códons carregada (%d entradas).", len(self._score))

    def analyze(self, gene: str, mrna: SeqRecord) -> dict:
        seq = str(mrna.seq).upper()
        cds = next(f for f in mrna.features if f.type == "CDS")
        start = int(cds.location.start)
        cds_seq = seq[start:int(cds.location.end)]
        kozak   = seq[max(0, start - 5): start + 4]  # -5…+3 (9 nt)

        scores = [self._score[c] for c in re.findall(r"...", cds_seq) if c in self._score]
        if not scores:
            raise ValueError("Nenhum códon com pontuação")

        gm = round(geometric_mean(scores), 4)
        return {
            "gene": gene,
            "cds_start": start,
            "kozak_seq": kozak,
            "geom_mean_cds": gm,
        }
