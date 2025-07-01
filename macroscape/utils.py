import re
import time
import logging
from contextlib import contextmanager

log = logging.getLogger("utils")


def snake_case(s: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^0-9a-zA-Z]+", "_", s)).strip("_").lower()


@contextmanager
def timed(msg: str):
    log.info(msg + "…")
    t0 = time.perf_counter()
    try:
        yield
    finally:
        log.info("%s em %.1fs", msg, time.perf_counter() - t0)
