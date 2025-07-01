import logging
import socket
import time
from http.client import HTTPException
from urllib.error import HTTPError, URLError
from typing import Protocol

from Bio import Entrez, SeqIO
from tenacity import retry, wait_exponential, stop_after_attempt, retry_if_exception

from .. import settings  # singleton
from ..utils import timed

log = logging.getLogger("ncbi")
socket.setdefaulttimeout(settings.SOCKET_TIMEOUT)

Entrez.email = settings.ENTREZ_EMAIL
Entrez.api_key = settings.ENTREZ_API_KEY
_NCBI_DELAY = settings.NCBI_DELAY


def _retryable(exc: Exception) -> bool:
    from requests.exceptions import RequestException
    if isinstance(exc, HTTPError):
        log.warning("HTTP %s – %s", getattr(exc, "code", "?"),
                                   getattr(exc, "reason", exc))
    return isinstance(exc, (HTTPError, URLError, RequestException,
                            socket.timeout, TimeoutError, HTTPException))


class INcbiClient(Protocol):
    """Interface (útil p/ testes)."""

    def get_mrna(self, gene: str): ...

    def get_protein_seq(self, gene: str) -> str: ...


class EntrezClient(INcbiClient):
    """Wrapper resiliente para E‑utilities."""

    @retry(stop=stop_after_attempt(5),
           wait=wait_exponential(multiplier=1, min=2, max=32),
           retry=retry_if_exception(_retryable))
    def _esearch(self, term: str, db: str = "nucleotide") -> str | None:
        h = Entrez.esearch(db=db, term=term, retmax=1, retmode="xml")
        rec = Entrez.read(h);
        h.close()
        time.sleep(_NCBI_DELAY)
        return rec.get("IdList", [None])[0]

    def _fetch_gb(self, uid: str, db: str):
        """Tenta gb → gbwithparts.  Levanta LookupError se nada disponível."""
        last_exc = None
        for rettype in ("gb", "gbwithparts"):
            try:
                h = Entrez.efetch(db=db, id=uid, rettype=rettype, retmode="text")
                rec = SeqIO.read(h, "genbank")
                h.close();
                time.sleep(_NCBI_DELAY)
                return rec
            except Exception as e:
                last_exc = e
                log.warning("efetch %s falhou (%s)", rettype, e)
                try:
                    h.close()
                except Exception:
                    pass
        raise LookupError(f"GenBank indisponível para UID ({last_exc})")

    # ------------------ API pública ------------------

    def get_mrna(self, gene: str):
        # 1ª estratégia: transcript variant 1
        term1 = f"Mus musculus[ORGN] AND {gene}[Gene] AND transcript variant 1[Title] AND mRNA[Filter]"
        uid = self._esearch(term1)
        if not uid:
            uid = self._esearch(f"Mus musculus[ORGN] AND {gene}[Gene] AND mRNA[Filter]")
        if not uid:
            raise LookupError("UID não encontrado")
        return self._fetch_gb(uid, "nucleotide")

    def get_protein_seq(self, gene: str) -> str:
        uid = self._esearch(f"Mus musculus[ORGN] AND {gene}[Gene]", db="protein")
        if not uid:
            raise LookupError("UID não encontrado")
        # UID → accession
        h = Entrez.esummary(db="protein", id=uid, retmode="xml")
        summ = Entrez.read(h);
        h.close()
        acc = summ[0].get("AccessionVersion") or summ[0].get("Caption")
        time.sleep(_NCBI_DELAY)
        # accession → sequência
        h = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
        seq = SeqIO.read(h, "fasta").seq;
        h.close()
        time.sleep(_NCBI_DELAY)
        return str(seq)
