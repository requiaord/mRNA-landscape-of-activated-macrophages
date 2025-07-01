# macroscape/main.py
import argparse
from pathlib import Path

from . import settings
from .pipeline import run_kozak, run_charge


def _add_common_args(p):
    p.add_argument("--excel",       type=Path, help="Planilha Excel")
    p.add_argument("--codon-table", type=Path, help="CSV de códons")
    p.add_argument("--win",         type=int,  help="Janela (aa) p/ análise de carga")


def main():
    parser = argparse.ArgumentParser(
        description="Macroscape – Kozak & Charge analyses"
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_all = sub.add_parser("all", help="Executa Kozak + Cargas")
    _add_common_args(p_all)

    p_kozak = sub.add_parser("kozak", help="Somente análise de Kozak")
    _add_common_args(p_kozak)

    p_charge = sub.add_parser("charge", help="Somente análise de carga")
    _add_common_args(p_charge)

    args = parser.parse_args()

    # valores ⇢ fallback para Settings
    excel   = args.excel       or settings.EXCEL_FILE
    table   = args.codon_table or settings.CODON_TABLE
    win     = args.win         or settings.WINDOW_SIZE

    if args.cmd in ("all", "kozak"):
        run_kozak(excel, table)
    if args.cmd in ("all", "charge"):
        run_charge(excel, win)


if __name__ == "__main__":
    main()
