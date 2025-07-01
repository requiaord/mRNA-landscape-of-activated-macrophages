from pathlib import Path
import pandas as pd
from ..utils import snake_case


class ExcelGeneLoader:
    """Carrega listas de genes (uma coluna = uma lista)."""

    def __init__(self, path: Path, sheet: int = 0):
        self._path = path
        self._sheet = sheet

    def load(self) -> dict[str, list[str]]:
        df = pd.read_excel(self._path, sheet_name=self._sheet)
        df.columns = map(snake_case, df.columns)

        for col in df.select_dtypes("object"):
            df[col] = (
                df[col]
                .fillna("")
                .astype(str)
                .str.split(r"[;,]")
                .apply(lambda lst: [g.strip() for g in lst if g.strip()])
            )
        return {c: [g for sub in df[c] for g in sub]
                for c in df.select_dtypes("object")}
