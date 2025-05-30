SEX_MAP = {"female": 0, "male": 1, "unknown": -1}

def sex_to_numeric(label: str) -> int:
    """'female'/'male'/'unknown' → 0/1/-1."""
    return SEX_MAP[label.lower()]

def numeric_to_sex(code: int) -> str:
    """0/1/-1 → 'female'/'male'/'unknown'."""
    inv = {v: k for k, v in SEX_MAP.items()}
    return inv[code]
