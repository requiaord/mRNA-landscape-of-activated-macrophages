import logging
from datetime import datetime
from pathlib import Path

from . import settings
from .io.excel_loader import ExcelGeneLoader
from .io.csv_writer import write_csv
from .ncbi.entrez_client import EntrezClient
from .analysis.codon import KozakAnalyzer
from .analysis.charge import ChargeAnalyzer
from .utils import timed, snake_case

log = logging.getLogger("pipeline")
_TS = datetime.now().strftime("%Y%m%d_%H%M%S")


def run_kozak(excel: Path, codon_table: Path):
    client   = EntrezClient()
    analyzer = KozakAnalyzer(codon_table)
    lists    = ExcelGeneLoader(excel).load()

    for col, genes in lists.items():
        ok, err = [], []
        log.info("→ Kozak: lista %s (%d genes)", col, len(genes))
        for g in genes:
            try:
                rec = client.get_mrna(g)
                ok.append(analyzer.analyze(g, rec))
            except Exception as e:
                err.append({"gene": g, "error": str(e)})
                log.warning("%s → %s", g, e)
        base = snake_case(col)
        write_csv(settings.PRODUCT_DIR / f"{base}_kozak_{_TS}.csv",
                  ["gene", "cds_start", "geom_mean_cds", "kozak_seq"],
                  [list(r.values()) for r in ok])
        write_csv(settings.PRODUCT_DIR / f"{base}_kozak_err_{_TS}.csv",
                  ["gene", "error"],
                  [list(r.values()) for r in err])


def run_charge(excel: Path, win: int):
    client   = EntrezClient()
    analyzer = ChargeAnalyzer()
    lists    = ExcelGeneLoader(excel).load()

    for col, genes in lists.items():
        rows30, rows100, err = [], [], []
        log.info("→ Charge: lista %s (%d genes)", col, len(genes))
        for g in genes:
            try:
                seq = client.get_protein_seq(g)
                res = analyzer.analyze(g, seq, win)
                rows30.append([g, res["seq_length"], *res["charges"]])
                rows100.append([g, res["seq_length"], *res["windows"]])
            except Exception as e:
                err.append([g, str(e)])
                log.warning("%s → %s", g, e)
        base = snake_case(col)
        hdr30 = ["gene", "seq_length"] + [f"R{i+1}" for i in range(win)]
        hdr100 = ["gene", "seq_length"] + [f"W{i+1}" for i in range(100)]
        write_csv(settings.PRODUCT_DIR / f"{base}_{win}res_{_TS}.csv", hdr30, rows30)
        write_csv(settings.PRODUCT_DIR / f"{base}_{win}win_{_TS}.csv", hdr100, rows100)
        write_csv(settings.PRODUCT_DIR / f"{base}_charge_err_{_TS}.csv",
                  ["gene", "error"], err)
