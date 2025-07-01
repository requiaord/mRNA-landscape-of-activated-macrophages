"""Convenience re‑exports."""
from .settings import Settings
from .log import setup as _setup_logging

settings = Settings()
console  = _setup_logging()
