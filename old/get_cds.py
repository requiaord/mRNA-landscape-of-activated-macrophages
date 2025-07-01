#!/usr/bin/env python3
import logging
from Bio import Entrez, SeqIO

# ------------------------------------------------------------------
# Configuração
# ------------------------------------------------------------------
Entrez.email = "seu_email@exemplo.com"  # ← substitua pelo seu e‑mail
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger("fetch_il6_variant1")

def fetch_il6_variant1_record() -> SeqIO.SeqRecord:
    """
    Busca e retorna o registro GenBank de mRNA para o gene Il6,
    transcript variant 1, em Mus musculus.
    """
    search_term = (
        "Mus musculus [ORGN] AND Il6 [Gene] "
        "AND transcript variant 1 [Title] AND mRNA [Filter]"
    )
    logger.info(f"Searching NCBI for variant 1 with: {search_term}")
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmode="xml", retmax=1)
    rec = Entrez.read(handle); handle.close()
    ids = rec.get("IdList", [])
    if not ids:
        raise ValueError("Nenhum ID encontrado para Il6 transcript variant 1")
    seq_id = ids[0]
    logger.info(f"Found ID for variant 1: {seq_id}")

    handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
    gb = SeqIO.read(handle, "genbank"); handle.close()
    return gb

def extract_cds_info(record: SeqIO.SeqRecord):
    """
    Encontra a primeira feature 'CDS' e retorna (start, end, length).
    """
    for f in record.features:
        if f.type == "CDS":
            start = int(f.location.start)
            end   = int(f.location.end)
            length = end - start
            logger.info(f"CDS located: start={start}, end={end}, length={length}")
            return start, end, length
    raise ValueError("Feature CDS não encontrada")

def main():
    try:
        rec = fetch_il6_variant1_record()
        full_seq = str(rec.seq).upper()
        logger.info(f"Full mRNA length: {len(full_seq)} nt")

        cds_start, cds_end, cds_len = extract_cds_info(rec)

        # 5' UTR
        utr5_len = cds_start
        utr5_seq = full_seq[:cds_start]

        # Região Kozak (5 nt upstream + AUG + 1 nt downstream = 9 nt)
        k_start = max(0, cds_start - 5)
        k_end   = cds_start + 4
        kozak_seq = full_seq[k_start:k_end]
        logger.info(f"Kozak region slice: {k_start}:{k_end} (length {len(kozak_seq)})")

        # Tradução do CDS
        cds_seq = rec.seq[cds_start:cds_end]
        aa_seq  = cds_seq.translate(to_stop=True)

        # Impressão dos resultados
        print("=== IL6 (Mus musculus) mRNA variant 1 ===")
        print(f"Accession         : {rec.id}")
        print(f"Description       : {rec.description}")
        print(f"mRNA total length : {len(full_seq)} nt")
        print(f"CDS start pos     : {cds_start}")
        print(f"CDS end pos       : {cds_end}")
        print(f"CDS length        : {cds_len} nt")
        print(f"\n5' UTR length     : {utr5_len} nt")
        print(f"5' UTR sequence   : {utr5_seq}")
        print(f"\nKozak region ({k_start}:{k_end}) length {len(kozak_seq)}: {kozak_seq}")
        print("\nTranslated CDS (AA):")
        print(aa_seq)

    except Exception as e:
        logger.critical(f"Erro geral: {e}")

if __name__ == "__main__":
    main()
