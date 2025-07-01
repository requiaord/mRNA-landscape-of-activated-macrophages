def _charges(seq: str) -> list[int]:
    cmap = {'K': 1, 'R': 1, 'D': -1, 'E': -1}
    vals = [cmap.get(a, 0) for a in seq]
    if vals:
        vals[0] = 1
    if len(vals) > 1:
        vals[-1] = -1
    return vals


def _windows(vals: list[int], size: int):
    return [sum(vals[i:i + size]) for i in range(len(vals) - size + 1)]


class ChargeAnalyzer:
    """Cargas residuais e somas em janelas deslizantes."""

    def analyze(self, gene: str, seq: str, win: int = 30, top: int = 100):
        ch = _charges(seq)
        win_sums = _windows(ch, win)
        return {
            "gene": gene,
            "seq_length": len(seq),
            "charges": ch[:win],        # primeiro N resíduos
            "windows": win_sums[:top],  # primeiras N janelas
        }
