#!/usr/bin/env python
"""
main.py

This script reads an Excel file containing multiple lists of genes (one list per column),
processes each list of genes by retrieving protein information from the NCBI (for Mus musculus),
computes per-residue charges and sliding window sums, and then outputs two CSV files for each gene list:
  1) A CSV with the first 30 individual residue charges (columns: Gene, SeqLength, R1..R30)
  2) A CSV with the first 100 sliding window sums of size 30 (columns: Gene, SeqLength, W1..W100)

Any genes that could not be processed (due to missing data, etc.) are logged in a separate error CSV.

Output CSV files are stored in the "product" folder and error logs in the "LOG" folder.
"""

import os
import re
import csv
import logging
import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from Bio import Entrez, SeqIO
from tenacity import retry, stop_after_attempt, wait_fixed, retry_if_exception_type

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
# Set your email (required by NCBI)
ENTREZ_EMAIL = "seu_email@exemplo.com"
Entrez.email = ENTREZ_EMAIL

# Input Excel file (assumed to be in the "data" folder)
EXCEL_FILE = os.path.join("../data", "Upregulated_proteins_BKP.xlsx")

# Output directories
PRODUCT_DIR = os.path.join(os.getcwd(), "product")
LOG_DIR = os.path.join(os.getcwd(), "LOG")
os.makedirs(PRODUCT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

# Timestamp for output filenames
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

# ------------------------------------------------------------------------------
# Logging Configuration
# ------------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("main")


# ------------------------------------------------------------------------------
# Helper functions for Excel column processing
# ------------------------------------------------------------------------------
def sanitize_column_name(col_name: str) -> str:
    """
    Convert the column name into a valid Python variable name by converting
    to lowercase and replacing non-alphanumeric characters with underscores.

    E.g., "All upregulated" -> "all_upregulated"
    """
    sanitized = re.sub(r'\W+', '_', col_name.lower())
    return sanitized


def process_column(series: pd.Series) -> list:
    """
    Processes a pandas Series by:
      - Converting values to strings.
      - Splitting values by ";" if present.
      - Stripping whitespace.
      - Removing empty items.
      - Removing duplicates.

    Returns:
        Sorted list of unique items.
    """
    items_set = set()
    for val in series.dropna().astype(str):
        parts = val.split(";") if ";" in val else [val]
        for part in parts:
            item = part.strip()
            if item:
                items_set.add(item)
    return sorted(items_set)


# ------------------------------------------------------------------------------
# Protein processing functions with retry decorators
# ------------------------------------------------------------------------------
@retry(stop=stop_after_attempt(3), wait=wait_fixed(2),
       retry=retry_if_exception_type(requests.exceptions.RequestException))
def search_protein_for_gene(gene_name: str, organism: str = "mus musculus") -> str:
    """
    Searches the NCBI Protein database for a given gene.
    Constructs the search term "organism [ORGN] AND gene_name" and returns the first UID found.

    Parameters:
        gene_name (str): Gene name.
        organism (str): Organism (default: "mus musculus").

    Returns:
        str: The first UID found, or None.
    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    term = f"{organism} [ORGN] AND {gene_name}"
    params = {"db": "protein", "term": term, "retmode": "json", "retmax": 1}
    logger.info(f"Searching UID for gene '{gene_name}' using term '{term}'")
    response = requests.get(url, params=params)
    response.raise_for_status()
    data = response.json()
    uid_list = data.get("esearchresult", {}).get("idlist", [])
    if not uid_list:
        logger.warning(f"No UID found for gene '{gene_name}'")
        return None
    uid = uid_list[0]
    logger.info(f"Found UID '{uid}' for gene '{gene_name}'")
    return uid


@retry(stop=stop_after_attempt(3), wait=wait_fixed(2),
       retry=retry_if_exception_type(requests.exceptions.RequestException))
def get_accession_from_uid(uid: str) -> str:
    """
    Retrieves the accession version (e.g., NP_000537.1) for the given UID using esummary.

    Parameters:
        uid (str): NCBI UID.

    Returns:
        str: Accession version, or None if not found.
    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {"db": "protein", "id": uid, "retmode": "json"}
    logger.info(f"Retrieving accession for UID '{uid}'")
    response = requests.get(url, params=params)
    response.raise_for_status()
    data = response.json()
    summary = data.get("result", {}).get(uid, {})
    accession = summary.get("accessionversion", None)
    if not accession:
        logger.warning(f"Accession not found for UID '{uid}'")
    return accession


@retry(stop=stop_after_attempt(3), wait=wait_fixed(2), retry=retry_if_exception_type(Exception))
def fetch_sequence(accession: str) -> tuple[str, str]:
    """
    Fetches the protein sequence and description for a given accession from NCBI.

    Parameters:
        accession (str): Protein accession.

    Returns:
        tuple[str, str]: (Sequence string, description)
    """
    logger.info(f"Fetching sequence for accession '{accession}'")
    handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    seq = str(record.seq)
    desc = record.description
    logger.info(f"Fetched sequence for '{accession}' (length: {len(seq)})")
    return seq, desc


def compute_charges(seq: str) -> list[int]:
    """
    Converts a protein sequence into charge values: +1 for K/R, -1 for D/E, 0 otherwise.
    Forces first residue to +1 and last to -1.

    Parameters:
        seq (str): Protein sequence.

    Returns:
        list[int]: List of charge values.
    """
    charge_map = {'K': 1, 'R': 1, 'D': -1, 'E': -1}
    charges = [charge_map.get(aa.upper(), 0) for aa in seq]
    if charges:
        charges[0] = 1
    if len(charges) > 1:
        charges[-1] = -1
    return charges


def sliding_window_sums(charges: list[int], window_size: int = 30) -> list[int]:
    """
    Computes sliding window sums of the charge values using the specified window size.

    Parameters:
        charges (list[int]): List of charge values.
        window_size (int): Window size (default: 30).

    Returns:
        list[int]: List of net charge sums for each window.
    """
    results = []
    n = len(charges)
    for start in range(n - window_size + 1):
        results.append(sum(charges[start:start + window_size]))
    return results


# ------------------------------------------------------------------------------
# ProteinProcessor Class
# ------------------------------------------------------------------------------
class ProteinProcessor:
    """
    Encapsulates the retrieval and processing of protein data for a given gene.
    """

    def __init__(self, gene: str, organism: str = "mus musculus"):
        """
        Initialize the processor with a gene name and organism.

        Parameters:
            gene (str): The gene name.
            organism (str): The organism (default "mus musculus").
        """
        self.gene = gene
        self.organism = organism
        self.uid = None
        self.accession = None
        self.sequence = None
        self.description = None
        self.seq_length = None
        self.charges = None
        self.windows = None

    def process(self) -> None:
        """
        Executes the full processing pipeline for the gene:
          - Searches for UID.
          - Retrieves the accession.
          - Fetches the sequence and description.
          - Computes individual residue charges.
          - Computes sliding window sums (window size 30).
        """
        self.uid = search_protein_for_gene(self.gene, self.organism)
        if not self.uid:
            logger.warning(f"No UID found for gene '{self.gene}'")
            return
        self.accession = get_accession_from_uid(self.uid)
        if not self.accession:
            logger.warning(f"No accession found for gene '{self.gene}' (UID: {self.uid})")
            return
        self.sequence, self.description = fetch_sequence(self.accession)
        self.seq_length = len(self.sequence)
        self.charges = compute_charges(self.sequence)
        self.windows = sliding_window_sums(self.charges, 30)

    def get_100_windows_row(self) -> list:
        """
        Returns a CSV row for the first 100 sliding window sums.

        Format: [Gene, SeqLength, W1, W2, ..., W100]
        If fewer than 100 windows exist, remaining cells are left empty.
        """
        row = [self.gene, self.seq_length if self.seq_length is not None else ""]
        for i in range(100):
            if self.windows and i < len(self.windows):
                row.append(self.windows[i])
            else:
                row.append("")
        return row

    def get_30_residues_row(self) -> list:
        """
        Returns a CSV row for the first 30 individual residue charge values.

        Format: [Gene, SeqLength, R1, R2, ..., R30]
        If the sequence has fewer than 30 residues, the remaining cells are left empty.
        """
        row = [self.gene, self.seq_length if self.seq_length is not None else ""]
        for i in range(30):
            if self.charges and i < len(self.charges):
                row.append(self.charges[i])
            else:
                row.append("")
        return row


# ------------------------------------------------------------------------------
# Function to process a gene (for parallel processing)
# ------------------------------------------------------------------------------
def process_gene(gene: str) -> tuple:
    """
    Processes a single gene and returns a tuple:
       (gene, row_100_windows, row_30_residues, error_message)
    If processing fails, error_message is not None.
    """
    processor = ProteinProcessor(gene)
    try:
        processor.process()
        if processor.sequence is None:
            return (gene, None, None, "No sequence found")
        return (gene, processor.get_100_windows_row(), processor.get_30_residues_row(), None)
    except Exception as e:
        return (gene, None, None, str(e))


# ------------------------------------------------------------------------------
# CSV writing helper
# ------------------------------------------------------------------------------
def write_csv(file_path: str, header: list, rows: list) -> None:
    """
    Writes a CSV file with the specified header and rows.

    Parameters:
        file_path (str): Path to the output CSV file.
        header (list): List of column headers.
        rows (list): List of rows (each row is a list).
    """
    try:
        with open(file_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(rows)
        logger.info(f"CSV file '{file_path}' generated successfully.")
    except Exception as e:
        logger.critical(f"Critical error writing '{file_path}': {e}")


# ------------------------------------------------------------------------------
# Main Function: Read Excel, Process Genes, and Write CSV Outputs
# ------------------------------------------------------------------------------
def main():
    # Read gene lists from Excel file. Assume each column is a separate list.
    excel_file = os.path.join("../data", "Upregulated_proteins_BKP.xlsx")
    try:
        df = pd.read_excel(excel_file)
        logger.info(f"Excel file '{excel_file}' loaded successfully.")
    except Exception as e:
        logger.critical(f"Failed to read Excel file '{excel_file}': {e}")
        return

    # For each column, process the gene list.
    # Use 'process_column' to get a unique gene list from each column.
    # We'll generate two CSVs per column.
    for col in df.columns:
        sanitized_col = sanitize_column_name(col)
        gene_list = process_column(df[col])
        logger.info(f"Processing column '{col}' as '{sanitized_col}' with {len(gene_list)} genes.")

        # Containers for CSV rows and errors for this gene list.
        results_100win = []  # For 100 sliding windows CSV.
        results_30res = []  # For 30 residues CSV.
        error_rows = []  # For errors: each row: [Gene, Error]

        # Process genes in parallel.
        with ThreadPoolExecutor(max_workers=8) as executor:
            future_to_gene = {executor.submit(process_gene, gene): gene for gene in gene_list}
            for future in as_completed(future_to_gene):
                gene = future_to_gene[future]
                try:
                    gene_name, row_100, row_30, err = future.result()
                    if err or row_100 is None or row_30 is None:
                        error_rows.append([gene_name, err if err else "Processing error"])
                    else:
                        results_100win.append(row_100)
                        results_30res.append(row_30)
                    logger.info(f"Processed gene: {gene_name}")
                except Exception as ex:
                    logger.error(f"Gene {gene} generated an exception: {ex}")
                    error_rows.append([gene, str(ex)])

        # Define output file paths for this column.
        csv_100_path = os.path.join(PRODUCT_DIR, f"{sanitized_col}_100windows_{TIMESTAMP}.csv")
        csv_30_path = os.path.join(PRODUCT_DIR, f"{sanitized_col}_30residues_{TIMESTAMP}.csv")
        error_csv = os.path.join(LOG_DIR, f"{sanitized_col}_errors_{TIMESTAMP}.csv")

        # Define headers.
        header_100win = ["Gene", "SeqLength"] + [f"W{i + 1}" for i in range(100)]
        header_30res = ["Gene", "SeqLength"] + [f"R{i + 1}" for i in range(30)]
        error_header = ["Gene", "Error"]

        # Write CSV files.
        write_csv(csv_100_path, header_100win, results_100win)
        write_csv(csv_30_path, header_30res, results_30res)
        write_csv(error_csv, error_header, error_rows)


if __name__ == "__main__":
    from concurrent.futures import ThreadPoolExecutor, as_completed

    main()
