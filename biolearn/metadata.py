"""
Metadata utilities for Biolearn.

Standard encodings for sex / age and lightweight search of library.yaml.
"""

from __future__ import annotations


import math
from typing import Union


def standardize_sex(value: Union[str, int, float, None]) -> Union[int, float]:
    """
    Standardize sex values to the DNA Methylation Array Data Standard.

    Standard: 0 for female, 1 for male, NaN for unknown/missing values.

    Args:
        value: Input sex value in various formats

    Returns:
        0 for female, 1 for male, NaN for unknown/missing
    """
    if value is None:
        return float("nan")

    # Handle numeric inputs
    if isinstance(value, (int, float)):
        if math.isnan(value):
            return float("nan")
        # Convert common numeric encodings
        if value == 0:
            return 0  # female
        elif value == 1:
            return 1  # male
        elif value == 2:
            return 0  # female (GEO encoding: 1=male, 2=female)
        else:
            return float("nan")

    # Handle string inputs
    if isinstance(value, str):
        normalized = value.strip().lower()
        mapping = {
            "female": 0,
            "f": 0,
            "male": 1,
            "m": 1,
            "unknown": float("nan"),
            "nan": float("nan"),
            "": float("nan"),
        }
        return mapping.get(normalized, float("nan"))

    return float("nan")


def sex_to_numeric(label: str) -> Union[int, float]:
    """Convert sex label to numeric standard format."""
    return standardize_sex(label)


def numeric_to_sex(code: Union[int, float]) -> str:
    """Convert numeric sex code to string representation."""
    if isinstance(code, (int, float)) and math.isnan(code):
        return "unknown"
    elif code == 0:
        return "female"
    elif code == 1:
        return "male"
    else:
        return "unknown"


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
