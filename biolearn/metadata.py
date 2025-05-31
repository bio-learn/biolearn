"""
Metadata utilities for Biolearn.

Standard encodings for sex / age and lightweight search of library.yaml.
"""

from __future__ import annotations

# ---------------------------------------------------------------------
# 1. Sex / age standard helpers (Issue #122)
# ---------------------------------------------------------------------
SEX_MAP: dict[str, int] = {"female": 0, "male": 1, "unknown": -1}


def sex_to_numeric(label: str) -> int:
    """'female' → 0, 'male' → 1, 'unknown' → -1."""
    return SEX_MAP[label.lower()]


def numeric_to_sex(code: int) -> str:
    """0 → 'female', 1 → 'male', -1 → 'unknown'."""
    inv = {v: k for k, v in SEX_MAP.items()}
    return inv[code]


# ---------------------------------------------------------------------
# 2. Iterator that copes with every library.yaml layout
# ---------------------------------------------------------------------
def _iter_library_items():
    """
    Yield ``(series_id, entry_dict)`` pairs.

    Handles all known *library.yaml* formats, including mappings,
    datasets, and flat item arrays.
    """
    import os
    from importlib import resources
    from pathlib import Path

    import yaml

    lib_path = Path(
        os.getenv("BIOLEARN_LIBRARY_PATH")
        or resources.files("biolearn.data").joinpath("library.yaml")
    )
    lib = yaml.safe_load(lib_path.read_text())

    # layout A: mapping at top level
    if isinstance(lib, dict) and "metadata" in next(iter(lib.values())):
        yield from lib.items()
        return

    # layout B: wrapped mapping
    if isinstance(lib, dict) and "datasets" in lib:
        yield from lib["datasets"].items()
        return

    # layout C: array of items
    if isinstance(lib, dict) and "items" in lib:
        for entry in lib["items"]:
            sid = entry.get("id") or entry.get("series_id")
            if sid:
                yield sid, entry
        return

    raise ValueError("Unrecognised library.yaml layout")


# ---------------------------------------------------------------------
# 3. Public helper: search without downloading big matrices (Issue #44)
# ---------------------------------------------------------------------
def search_metadata(**criteria):
    """
    Return a filtered DataFrame based on metadata criteria.

    Examples:
        >>> search_metadata(sex="male", min_age=70)
        >>> search_metadata(sex="female")
    """
    import pandas as pd

    hits: list[dict] = []

    for sid, entry in _iter_library_items():
        # -------- resolve a flat metadata dict ------------------------
        if "metadata" in entry:
            meta = entry["metadata"]
        else:
            meta: dict = {}
            parser = entry.get("parser", {})

            # a) parser-block values
            if "sex" in parser and "parse" in parser["sex"]:
                meta["sex"] = parser["sex"]["parse"]
            if "age" in parser and "parse" in parser["age"]:
                meta["age"] = float(parser["age"]["parse"])

            # b) direct top-level fallbacks
            if "sex" in entry:
                meta.setdefault("sex", entry["sex"])
            if "age" in entry and isinstance(entry["age"], (int, float, str)):
                meta.setdefault("age", float(entry["age"]))
        # ----------------------------------------------------------------

        # -------- sex filter (robust normalisation) --------------------
        if "sex" in criteria:
            wanted = criteria["sex"].lower()
            raw = meta.get("sex")

            # numeric codes used by GEO/dbGaP
            if isinstance(raw, (int, float)):
                raw = {
                    0: "unknown",
                    1: "male",
                    2: "female",
                    -1: "unknown",
                }.get(int(raw))

            # strings / single letters
            if isinstance(raw, str):
                raw = {"f": "female", "m": "male"}.get(
                    raw.strip().lower(), raw.strip().lower()
                )

            # only skip if we actually know the sex and it mismatches
            if raw is not None and raw != wanted:
                continue
        # ----------------------------------------------------------------

        # -------- age filter -------------------------------------------
        if (
            "min_age" in criteria
            and float(meta.get("age", -1)) < criteria["min_age"]
        ):
            continue

        hits.append({"series_id": sid, **meta})

    return pd.DataFrame(hits)
