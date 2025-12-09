"""Post-load correction functions for GEO datasets with known issues."""

import pandas as pd

from biolearn.util import cached_download


def fix_gse110554(geo_data, source_path):
    """Fix cell_type metadata for GSM2998097 and GSM2998106.

    These samples have metadata on row 14 instead of row 13 due to
    GEO formatting issues. See https://github.com/bio-learn/biolearn/issues/87
    """
    samples_to_fix = ["GSM2998097", "GSM2998106"]
    file_path = cached_download(source_path)

    # Read correct cell_type values from row 53 (1-indexed: row 54)
    # The standard cell_type is on row 53, but these samples need row 54
    raw = pd.read_table(file_path, index_col=0, skiprows=52, nrows=1)

    for sample in samples_to_fix:
        if sample in raw.columns and sample in geo_data.metadata.index:
            value = raw[sample].iloc[0]
            if isinstance(value, str) and ":" in value:
                value = value.split(":")[1].strip()
            geo_data.metadata.loc[sample, "cell_type"] = value

    return geo_data


# Registry maps correction names to functions
CORRECTIONS = {
    "fix_gse110554": fix_gse110554,
}


def apply_correction(name, geo_data, source_path):
    """Apply a named correction to GeoData."""
    if name not in CORRECTIONS:
        raise ValueError(f"Unknown correction: {name}")
    return CORRECTIONS[name](geo_data, source_path)
