"""Unit conversion utilities for clinical biomarker data."""

import warnings
import pandas as pd
from biolearn.clinical.registry import BIOMARKER_REGISTRY


def convert_units(df, source_units=None, units=None):
    """Convert biomarker columns to canonical units.

    Parameters
    ----------
    df : DataFrame
        DataFrame with biomarker columns (samples as rows).
    source_units : str, optional
        Named source preset (e.g. 'ukbiobank'). Applies preset unit
        mappings for all biomarkers from that source.
    units : dict, optional
        Per-biomarker unit overrides. Keys are biomarker names, values
        are source unit strings (e.g. ``{"creatinine": "umol/L"}``).
        Overrides any preset from ``source_units``.

    Returns
    -------
    DataFrame
        Copy of df with converted columns.

    Raises
    ------
    ValueError
        If a specified unit has no known conversion.
    """
    unit_map = {}
    if source_units is not None:
        unit_map.update(BIOMARKER_REGISTRY.get_source_preset(source_units))
    if units is not None:
        unit_map.update(units)

    if not unit_map:
        return df.copy()

    result = df.copy()
    for biomarker, source_unit in unit_map.items():
        if biomarker not in result.columns:
            continue
        if biomarker not in BIOMARKER_REGISTRY:
            warnings.warn(
                f"Biomarker '{biomarker}' not in registry, skipping conversion."
            )
            continue

        entry = BIOMARKER_REGISTRY.get(biomarker)
        if source_unit == entry["unit"]:
            continue  # already in canonical units

        if source_unit not in entry["conversions"]:
            raise ValueError(
                f"No conversion from '{source_unit}' to '{entry['unit']}' "
                f"for biomarker '{biomarker}'. "
                f"Known source units: {list(entry['conversions'].keys())}"
            )

        converter = entry["conversions"][source_unit]
        result[biomarker] = result[biomarker].apply(converter)

    return result


def validate_ranges(df, warn=True):
    """Check biomarker values against expected ranges.

    Parameters
    ----------
    df : DataFrame
        DataFrame with biomarker columns (samples as rows).
    warn : bool
        If True, emit warnings for out-of-range values.

    Returns
    -------
    dict
        Mapping of biomarker name to count of out-of-range values.
    """
    out_of_range = {}
    for col in df.columns:
        if col not in BIOMARKER_REGISTRY:
            continue
        lo, hi = BIOMARKER_REGISTRY.valid_range(col)
        mask = (df[col] < lo) | (df[col] > hi)
        count = mask.sum()
        if count > 0:
            out_of_range[col] = int(count)
            if warn:
                warnings.warn(
                    f"Biomarker '{col}': {count} values outside expected "
                    f"range [{lo}, {hi}] (unit: {BIOMARKER_REGISTRY.canonical_unit(col)}). "
                    f"Check units."
                )
    return out_of_range
