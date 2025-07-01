import csv
import logging
from pathlib import Path

log = logging.getLogger("csv_writer")


def write_csv(path: Path, header: list[str], rows: list):
    path.parent.mkdir(exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    log.info("CSV salvo → %s (%d linhas)", path, len(rows))
