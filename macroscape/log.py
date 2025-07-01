import logging
import logging.handlers
from pathlib import Path
from rich.console import Console

_LOG_FILE = Path(__file__).resolve().parent / "logs" / "run.log"
_LOG_FILE.parent.mkdir(exist_ok=True)


def setup(level: str = "INFO") -> Console:
    console = Console()
    lvl     = getattr(logging, level.upper(), logging.INFO)

    # Handler colorido (stdout)
    sh = logging.StreamHandler()
    sh.setLevel(lvl)
    sh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s",
                                      "%H:%M:%S"))

    # Handler de arquivo rotativo
    fh = logging.handlers.RotatingFileHandler(
        _LOG_FILE, maxBytes=2_000_000, backupCount=3, encoding="utf-8"
    )
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
        "%Y-%m-%d %H:%M:%S"))

    logging.basicConfig(level=lvl, handlers=[sh, fh], force=True)
    console.log(f"[bold green]Log level → [yellow]{logging.getLevelName(lvl)}[/]")
    return console
