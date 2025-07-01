#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import re
from pathlib import Path


def snake_case(name: str) -> str:
    """
    Converte texto livre para snake_case
    """
    name = re.sub(r'[^0-9a-zA-Z]+', '_', name)
    name = re.sub(r'_{2,}', '_', name)
    return name.strip('_').lower()


def load_and_clean(filepath, sheet_name=0):
    """
    Lê e limpa o Excel de proteínas upregulated:
      - lê a aba especificada (0 -> primeira aba)
      - renomeia colunas para snake_case
      - split em ',' ou ';' em cada célula de texto
    """
    df = pd.read_excel(filepath, sheet_name=sheet_name)
    df.columns = [snake_case(col) for col in df.columns]

    for col in df.select_dtypes(include='object'):
        df[col] = (
            df[col]
            .fillna('')
            .astype(str)
            .str.split(r'[;,]')
            .apply(lambda lst: [g.strip() for g in lst if g.strip()])
        )
    return df


def flatten_columns(df):
    """
    Para cada coluna com listas, retorna um dicionário de listas achatadas.
    """
    flat = {}
    for col in df.select_dtypes(include='object'):
        # concatena todas as sublistas
        flat[col] = [gene for sub in df[col] for gene in sub]
    return flat


if __name__ == '__main__':
    default_file = Path(__file__).parent / 'data' / 'upregulated_proteins.xlsx'
    df = load_and_clean(default_file)

    # Explode para testar (cada gene em nova linha)
    exploded = df.copy()
    for col in exploded.select_dtypes(include='object'):
        exploded = exploded.explode(col)
    print("DataFrame explodido (primeiras 5 linhas):")
    print(exploded.head(), "\n")

    # Lista achatada por coluna
    flat_lists = flatten_columns(df)
    for col, genes in flat_lists.items():
        print(f"{col} = {genes}\n")
