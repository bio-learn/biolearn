import os
import numpy as np
import pandas as pd
import pytest
from biolearn.data_library import GeoData


def _make_clinical_df():
    """Create a small clinical DataFrame (samples as rows)."""
    return pd.DataFrame(
        {
            "age": [45, 62, 38],
            "sex": [1, 0, 1],
            "albumin": [4.2, 3.8, 4.5],
            "creatinine": [0.9, 1.1, 0.8],
            "glucose": [5.1, 6.2, 4.8],
            "white_blood_cell_count": [6.5, 8.0, 5.5],
            "lymphocyte_percent": [30.0, 25.0, 35.0],
            "mean_cell_volume": [88.0, 92.0, 86.0],
            "red_blood_cell_distribution_width": [12.5, 14.0, 12.0],
            "alkaline_phosphate": [65.0, 80.0, 55.0],
        },
        index=["P1", "P2", "P3"],
    )


def test_geodata_with_clinical_layer():
    """GeoData accepts a clinical DataFrame."""
    clinical = pd.DataFrame(
        {"P1": [4.2, 0.9], "P2": [3.8, 1.1]},
        index=["albumin", "creatinine"],
    )
    metadata = pd.DataFrame({"age": [45, 62]}, index=["P1", "P2"])
    geo = GeoData(metadata=metadata, clinical=clinical)

    assert geo.clinical is not None
    assert list(geo.clinical.columns) == ["P1", "P2"]
    assert list(geo.clinical.index) == ["albumin", "creatinine"]
    assert geo.dnam is None


def test_geodata_clinical_defaults_to_none():
    """Clinical layer defaults to None for backward compatibility."""
    metadata = pd.DataFrame({"age": [45]}, index=["P1"])
    geo = GeoData(metadata=metadata)
    assert geo.clinical is None


def test_from_clinical_matrix_basic():
    """from_clinical_matrix separates metadata and transposes biomarkers."""
    df = _make_clinical_df()
    geo = GeoData.from_clinical_matrix(df)

    # Metadata should contain age and sex
    assert "age" in geo.metadata.columns
    assert "sex" in geo.metadata.columns
    assert len(geo.metadata) == 3

    # Clinical should be features-as-rows, samples-as-columns
    assert geo.clinical is not None
    assert set(geo.clinical.columns) == {"P1", "P2", "P3"}
    assert "albumin" in geo.clinical.index
    assert "creatinine" in geo.clinical.index

    # Metadata columns should NOT be in clinical
    assert "age" not in geo.clinical.index
    assert "sex" not in geo.clinical.index


def test_from_clinical_matrix_preserves_values():
    """Values survive the transpose correctly."""
    df = _make_clinical_df()
    geo = GeoData.from_clinical_matrix(df)

    assert geo.clinical.loc["albumin", "P1"] == 4.2
    assert geo.clinical.loc["creatinine", "P2"] == 1.1
    assert geo.metadata.loc["P1", "age"] == 45


def test_from_clinical_matrix_no_metadata_cols():
    """Works when input has no metadata columns."""
    df = pd.DataFrame(
        {"albumin": [4.2, 3.8], "creatinine": [0.9, 1.1]},
        index=["P1", "P2"],
    )
    geo = GeoData.from_clinical_matrix(df)

    assert len(geo.metadata.columns) == 0
    assert geo.clinical is not None
    assert "albumin" in geo.clinical.index


def test_from_clinical_matrix_unit_conversion():
    """Unit conversion via the units parameter works."""
    df = pd.DataFrame(
        {"creatinine": [79.56, 97.24]},  # umol/L values
        index=["P1", "P2"],
    )
    geo = GeoData.from_clinical_matrix(df, units={"creatinine": "umol/L"})

    # Should be converted to mg/dL (divide by 88.42)
    converted = geo.clinical.loc["creatinine", "P1"]
    assert abs(converted - 79.56 / 88.42) < 0.01


def test_from_clinical_matrix_source_preset():
    """source_units preset applies correct conversions."""
    df = pd.DataFrame(
        {
            "albumin": [42.0],  # g/L (UK Biobank)
            "creatinine": [88.42],  # umol/L
        },
        index=["P1"],
    )
    geo = GeoData.from_clinical_matrix(df, source_units="ukbiobank")

    # albumin: 42 g/L -> 4.2 g/dL
    assert abs(geo.clinical.loc["albumin", "P1"] - 4.2) < 0.01
    # creatinine: 88.42 umol/L -> 1.0 mg/dL
    assert abs(geo.clinical.loc["creatinine", "P1"] - 1.0) < 0.01


def test_copy_preserves_clinical():
    """GeoData.copy() deep-copies the clinical layer."""
    df = _make_clinical_df()
    geo = GeoData.from_clinical_matrix(df)
    geo_copy = geo.copy()

    # Modify original
    geo.clinical.iloc[0, 0] = -999

    # Copy should be unaffected
    assert geo_copy.clinical.iloc[0, 0] != -999


def test_save_load_roundtrip_with_clinical(tmp_path):
    """Clinical data survives save_csv / load_csv roundtrip."""
    df = _make_clinical_df()
    geo = GeoData.from_clinical_matrix(df)

    folder = str(tmp_path)
    geo.save_csv(folder, "test")

    # Verify clinical file was created
    assert os.path.exists(os.path.join(folder, "test_clinical.csv"))

    # Load it back
    loaded = GeoData.load_csv(folder, "test", validate=False)
    assert loaded.clinical is not None
    assert set(loaded.clinical.index) == set(geo.clinical.index)
    assert set(loaded.clinical.columns) == set(geo.clinical.columns)

    # Values should match
    pd.testing.assert_frame_equal(
        loaded.clinical.sort_index(axis=0).sort_index(axis=1),
        geo.clinical.sort_index(axis=0).sort_index(axis=1),
        atol=1e-10,
    )


def test_validate_metadata_omics_includes_clinical():
    """Validation includes clinical samples in consistency check."""
    clinical = pd.DataFrame(
        {"P1": [4.2], "P2": [3.8], "P3": [4.5]},
        index=["albumin"],
    )
    metadata = pd.DataFrame({"age": [45, 62]}, index=["P1", "P2"])
    geo = GeoData(metadata=metadata, clinical=clinical)

    with pytest.warns(UserWarning, match="without metadata"):
        geo._validate_metadata_omics_consistency()
