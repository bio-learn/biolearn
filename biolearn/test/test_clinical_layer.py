import os
import numpy as np
import pandas as pd
import pytest
from biolearn.data_library import GeoData


def _make_clinical_df():
    """Build a small DataFrame in the canonical samples-as-rows shape.

    Values are in biolearn canonical units (g/L for albumin, umol/L for
    creatinine, mmol/L for glucose, etc.).
    """
    return pd.DataFrame(
        {
            "age": [45, 62, 38],
            "sex": [1, 0, 1],
            "albumin": [42.0, 38.0, 45.0],
            "creatinine": [80.0, 97.0, 71.0],
            "glucose": [5.1, 6.2, 4.8],
            "white_blood_cell_count": [6.5, 8.0, 5.5],
            "lymphocyte_percent": [30.0, 25.0, 35.0],
            "mean_cell_volume": [88.0, 92.0, 86.0],
            "red_blood_cell_distribution_width": [12.5, 14.0, 12.0],
            "alkaline_phosphate": [65.0, 80.0, 55.0],
        },
        index=pd.Index(["P1", "P2", "P3"], name="id"),
    )


def test_geodata_accepts_clinical_layer():
    """GeoData stores clinical data with samples-as-rows."""
    clinical = pd.DataFrame(
        {"albumin": [4.2, 3.8], "creatinine": [0.9, 1.1]},
        index=["P1", "P2"],
    )
    metadata = pd.DataFrame({"age": [45, 62]}, index=["P1", "P2"])
    geo = GeoData(metadata=metadata, clinical=clinical)

    assert geo.clinical is not None
    assert list(geo.clinical.index) == ["P1", "P2"]
    assert list(geo.clinical.columns) == ["albumin", "creatinine"]
    assert geo.dnam is None


def test_clinical_defaults_to_none():
    """The clinical layer is opt-in and defaults to None."""
    metadata = pd.DataFrame({"age": [45]}, index=["P1"])
    geo = GeoData(metadata=metadata)
    assert geo.clinical is None


def test_from_clinical_matrix_splits_metadata_and_biomarkers():
    """from_clinical_matrix moves age/sex/mortality into metadata."""
    df = _make_clinical_df()
    geo = GeoData.from_clinical_matrix(df)

    # Metadata grabs the demographic and mortality columns when present
    assert "age" in geo.metadata.columns
    assert "sex" in geo.metadata.columns
    assert len(geo.metadata) == 3

    # Clinical keeps samples-as-rows
    assert geo.clinical is not None
    assert set(geo.clinical.index) == {"P1", "P2", "P3"}
    assert "albumin" in geo.clinical.columns
    assert "creatinine" in geo.clinical.columns

    # Metadata fields must not leak into clinical
    assert "age" not in geo.clinical.columns
    assert "sex" not in geo.clinical.columns


def test_from_clinical_matrix_preserves_values():
    """Values stay attached to the correct sample after the split."""
    df = _make_clinical_df()
    geo = GeoData.from_clinical_matrix(df)

    assert geo.clinical.loc["P1", "albumin"] == 42.0
    assert geo.clinical.loc["P2", "creatinine"] == 97.0
    assert geo.metadata.loc["P1", "age"] == 45


def test_from_clinical_matrix_without_metadata_columns():
    """Works when the input is biomarkers only."""
    df = pd.DataFrame(
        # Albumin in g/L, creatinine in umol/L (canonical units)
        {"albumin": [42.0, 38.0], "creatinine": [80.0, 97.0]},
        index=["P1", "P2"],
    )
    geo = GeoData.from_clinical_matrix(df)

    assert len(geo.metadata.columns) == 0
    assert geo.clinical is not None
    assert "albumin" in geo.clinical.columns


def test_from_clinical_matrix_per_biomarker_unit_override():
    """The ``units`` argument applies per-biomarker conversions."""
    # Source mg/dL glucose becomes mmol/L (canonical) using factor 0.05551
    df = pd.DataFrame(
        {"glucose": [90.0, 180.0]},
        index=["P1", "P2"],
    )
    geo = GeoData.from_clinical_matrix(df, units={"glucose": "mg/dL"})

    assert abs(geo.clinical.loc["P1", "glucose"] - 90.0 * 0.05551) < 1e-6
    assert abs(geo.clinical.loc["P2", "glucose"] - 180.0 * 0.05551) < 1e-6


def test_from_clinical_matrix_fhs_source_preset():
    """The ``fhs`` preset converts glucose from mg/dL to mmol/L."""
    df = pd.DataFrame(
        {"glucose": [90.0]},  # mg/dL in raw FHS Period 1
        index=["P1"],
    )
    geo = GeoData.from_clinical_matrix(df, source_units="fhs")

    # glucose: 90 mg/dL * 0.05551 ≈ 4.9959 mmol/L
    assert abs(geo.clinical.loc["P1", "glucose"] - 90.0 * 0.05551) < 1e-6


def test_copy_preserves_clinical():
    """GeoData.copy() deep-copies the clinical layer."""
    df = _make_clinical_df()
    geo = GeoData.from_clinical_matrix(df)
    geo_copy = geo.copy()

    # Mutate the original; the copy should not change
    geo.clinical.iloc[0, 0] = -999

    assert geo_copy.clinical.iloc[0, 0] != -999


def test_save_load_roundtrip_with_clinical(tmp_path):
    """Clinical data survives save_csv / load_csv unchanged."""
    df = _make_clinical_df()
    geo = GeoData.from_clinical_matrix(df)

    folder = str(tmp_path)
    geo.save_csv(folder, "test")

    assert os.path.exists(os.path.join(folder, "test_clinical.csv"))

    loaded = GeoData.load_csv(folder, "test", validate=False)
    assert loaded.clinical is not None
    assert set(loaded.clinical.index) == set(geo.clinical.index)
    assert set(loaded.clinical.columns) == set(geo.clinical.columns)

    pd.testing.assert_frame_equal(
        loaded.clinical.sort_index(axis=0).sort_index(axis=1),
        geo.clinical.sort_index(axis=0).sort_index(axis=1),
        atol=1e-10,
    )


def test_validate_metadata_omics_includes_clinical():
    """Consistency check picks up samples that live only in the clinical layer."""
    clinical = pd.DataFrame(
        {"albumin": [4.2, 3.8, 4.5]},
        index=["P1", "P2", "P3"],
    )
    metadata = pd.DataFrame({"age": [45, 62]}, index=["P1", "P2"])
    geo = GeoData(metadata=metadata, clinical=clinical)

    with pytest.warns(UserWarning, match="without metadata"):
        geo._validate_metadata_omics_consistency()
