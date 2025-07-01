from pathlib import Path
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """
    Camada única de configuração.
    Pode ser sobrescrita por .env ou variáveis de ambiente prefixadas com GT_*
    """
    model_config = SettingsConfigDict(env_file=".env", env_prefix="GT_")

    # ---- caminhos principais ----
    BASE_DIR: Path = Path(__file__).resolve().parent
    DATA_DIR: Path = BASE_DIR / "data"
    LOG_DIR: Path = BASE_DIR / "logs"
    PRODUCT_DIR: Path = BASE_DIR / "product"

    EXCEL_FILE: Path = DATA_DIR / "upregulated_proteins.xlsx"
    CODON_TABLE: Path = DATA_DIR / "tabela.csv"

    # ---- NCBI / rede ----
    ENTREZ_EMAIL: str = "seu_email@exemplo.com"
    ENTREZ_API_KEY: str = ""  # vazio = sem chave
    SOCKET_TIMEOUT: int = 30  # segundos
    NCBI_DELAY: float | None = None  # calculado depois

    # ---- parâmetros de análise ----
    WINDOW_SIZE: int = 30

    def __init__(self, **kw):
        super().__init__(**kw)
        # cria pastas se não existirem
        for p in (self.DATA_DIR, self.LOG_DIR, self.PRODUCT_DIR):
            p.mkdir(exist_ok=True)

        # delay adaptativo
        if self.NCBI_DELAY is None:
            self.NCBI_DELAY = 0.12 if self.ENTREZ_API_KEY else 0.4
