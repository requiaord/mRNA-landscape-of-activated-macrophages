
# Gene-Based Protein Net Charge Analysis and CSV Generation

## Overview

This project contains a Python script (`generate_csv.py`) designed to perform protein data retrieval and analysis based on gene lists extracted from an Excel file. The script is aimed at researchers in bioinformatics who wish to process lists of genes, retrieve corresponding protein sequences for _Mus musculus_ from the NCBI database, and compute detailed charge metrics. These metrics include:

-   **Individual Residue Charge Values**: The charge for each amino acid (residue) in a protein sequence.
    
-   **Sliding Window Net Charge**: The net charge computed over a sliding window of 30 consecutive amino acids across the sequence.
    

For each gene list (each column in the Excel file), the script outputs **two CSV files**:

1.  A CSV file containing the first 30 individual residue charge values.
    
2.  A CSV file containing the first 100 sliding window net charge sums.
    

Any genes that could not be processed (due to missing sequence data or errors) are logged in separate error CSV files within a designated log folder.

## File Structure

-   **data/Upregulated_proteins.xlsx**  
    The input Excel file containing multiple lists of gene names. Each column represents a separate gene list.
    
-   **generate_csv.py**  
    The main script that:
    
    -   Reads the Excel file.
        
    -   Processes each column to extract unique gene names.
        
    -   Uses parallel processing to retrieve and analyze protein data from the NCBI for each gene.
        
    -   Computes per-residue charges and sliding window net charge sums.
        
    -   Generates output CSV files for each gene list.
        
-   **product/**  
    Directory where the generated CSV files (containing analysis results) are saved.
    
-   **LOG/**  
    Directory where error log CSV files are saved. Each error log contains the genes that could not be processed along with corresponding error messages.
    
-   **requirements.txt**  
    (Not shown here) Should list project dependencies (e.g., `pandas`, `biopython`, `tenacity`, `requests`, `openpyxl`).
    

## Configuration

The script uses several configurable parameters defined at the top of the file:

-   **NCBI Email**:  
    Set your email using the `ENTREZ_EMAIL` variable. This is required by NCBI to track API usage.
    
-   **Excel File Location**:  
    The input Excel file is assumed to be located in the `data` folder as `"data/Upregulated_proteins.xlsx"`.
    
-   **Output Directories**:  
    CSV output files are saved to the `product` directory, and error logs are saved to the `LOG` directory. The script creates these directories if they do not exist.
    
-   **Timestamp**:  
    A timestamp is generated to append to output file names, ensuring uniqueness.
    

## Data Processing Workflow

### 1. Reading the Excel File and Extracting Gene Lists

-   **Excel Reading**:  
    The script reads the provided Excel file using `pandas.read_excel()`. It assumes that each column represents a separate gene list.
    
-   **Column Processing**:  
    For each column, the following occurs:
    
    -   **Sanitization**: The column name is converted to a valid Python variable name (lowercase with underscores).
        
    -   **Value Processing**: Each cell value is converted to a string. If the string contains semicolons (`;`), it is split into multiple items. Whitespace is stripped, empty items are removed, and duplicate values are eliminated.
        
    -   **Unique Gene List**: A sorted list of unique gene names is created from each column.
        

### 2. Protein Data Retrieval and Analysis

For each gene in a list, the script processes the gene data in parallel using Python’s `ThreadPoolExecutor`:

-   **NCBI Search**:
    
    -   The script uses the NCBI E-utilities (specifically, the `esearch` endpoint) to search the Protein database with the query:  
        `mus musculus [ORGN] AND <gene_name>`.
        
    -   The first UID returned is used for further processing.
        
-   **Accession Retrieval**:
    
    -   Using the UID obtained from the search, the script calls the `esummary` endpoint to retrieve the accession version (e.g., `NP_000537.1`).
        
-   **Sequence Fetching**:
    
    -   The protein sequence and description are fetched using the `efetch` endpoint from NCBI.
        
-   **Charge Calculation**:
    
    -   **Per-Residue Charges**:  
        Each amino acid in the protein sequence is assigned a charge value:
        
        -   **+1** for lysine (K) and arginine (R).
            
        -   **-1** for aspartate (D) and glutamate (E).
            
        -   **0** for all other amino acids. The script enforces that the first residue is +1 and the last is -1.
            
    -   **Sliding Window Net Charge**:  
        The script calculates the net charge using a sliding window of 30 residues. For each window (e.g., positions 1–30, 2–31, etc.), the charge values of the residues within the window are summed. Only the first 100 windows are retained for the CSV output.
        

### 3. CSV File Generation

For each gene list (i.e., each column in the Excel file), the script generates two CSV files:

-   **CSV File 1: 100 Windows**
    
    -   **Columns**: `Gene, SeqLength, W1, W2, ..., W100`
        
    -   **Data**:
        
        -   **Gene**: The gene name.
            
        -   **SeqLength**: The total length (number of amino acids) of the protein sequence.
            
        -   **W1–W100**: The net charge sums for the first 100 sliding windows of size 30. If a protein yields fewer than 100 windows, the remaining entries are left blank.
            
-   **CSV File 2: 30 Residues**
    
    -   **Columns**: `Gene, SeqLength, R1, R2, ..., R30`
        
    -   **Data**:
        
        -   **Gene**: The gene name.
            
        -   **SeqLength**: The total length of the protein sequence.
            
        -   **R1–R30**: The individual charge values of the first 30 residues of the protein sequence. If a protein has fewer than 30 residues, the remaining entries are left blank.
            
-   **Error Logging**:  
    If any gene cannot be processed (e.g., no sequence found, retrieval error), the gene and its error message are logged in a separate CSV file in the `LOG` directory. The error log file is named based on the corresponding gene list’s column name and timestamp.
    

### 4. Parallel Processing

To speed up the processing of potentially large gene lists, the script uses a `ThreadPoolExecutor` (with a configurable number of workers, here set to 8) to process genes in parallel. Each gene is processed independently, and errors in processing one gene do not stop the overall pipeline.

## How to Run

1.  **Ensure Directory Structure**:  
    Place the Excel file (`Upregulated_proteins.xlsx`) in the `data` folder.
    
2.  **Install Dependencies**:  
    Make sure you have the required packages installed. For example:
    
    bash
    
    Copiar
    
    `pip install pandas openpyxl biopython tenacity requests` 
    
3.  **Execute the Script**:  
    Run the script from the project root:
    
    bash
    
    Copiar
    
    `python generate_csv.py` 
    
    The output CSV files will be saved in the `product` folder, and error logs (if any) in the `LOG` folder.
    

## File Naming

The output file names incorporate the sanitized column name from the Excel file and a timestamp. For example, you might see:

-   `upregulated_100windows_20220415_153210.csv`
    
-   `upregulated_30residues_20220415_153210.csv`
    

## Summary

-   The script reads multiple gene lists from an Excel file.
    
-   For each gene, it retrieves protein information from NCBI and computes:
    
    -   Per-residue charge values (with special handling of the first and last residues).
        
    -   Sliding window sums (net charge over windows of 30 residues), retaining the first 100 windows.
        
-   Two CSV files per gene list are generated:
    
    -   One with the first 100 sliding window net charge values.
        
    -   One with the first 30 individual residue charge values.
        
-   Errors during processing are logged into separate CSV files in the `LOG` folder.
    
-   The processing is done in parallel to increase efficiency.
    

This design employs modular, clean, and scalable code that adheres to SOLID principles, making it easier to maintain and extend for future bioinformatics analyses.