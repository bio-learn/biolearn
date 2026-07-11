from biolearn import load
from biolearn.data_library import DataLibrary, GeoData
import pytest

load_columns = ["sex", "age", "glucose", "is_dead", "months_until_death"]


def test_fhs_columns():
    df = load.load_fhs()
    verify_expected_columns(df)


def test_nhanes_columns():
    df = load.load_nhanes(2010)
    verify_expected_columns(df)


def test_can_load_nhanes_2012():
    df = load.load_nhanes(2012)


def test_expected_error_when_loading_unsupported_year_nhanes():
    with pytest.raises(ValueError):
        df = load.load_nhanes(1913)


def test_nhanes_2010_via_datalibrary_returns_clinical_layer():
    geo = DataLibrary().get("NHANES_2010").load()
    assert isinstance(geo, GeoData)
    assert geo.clinical is not None
    # Samples-as-rows orientation: each row is a participant
    assert "age" in geo.metadata.columns
    # Biomarkers live as columns on the clinical frame
    assert "glucose" in geo.clinical.columns


def test_fhs_via_datalibrary_applies_fhs_unit_conversion():
    """End-to-end check that the FHS source preset produces canonical units."""
    geo = DataLibrary().get("FHS").load()
    assert isinstance(geo, GeoData)
    assert geo.clinical is not None
    assert "glucose" in geo.clinical.columns

    # Raw FHS glucose is mg/dL (~70-200). After the fhs preset converts to
    # canonical mmol/L the median should be in the typical 4-12 range.
    glucose = geo.clinical["glucose"].dropna()
    assert (
        glucose.median() < 20
    ), "Glucose median above 20 suggests conversion didn't run (still mg/dL?)"
    assert (
        glucose.median() > 2
    ), "Glucose median below 2 suggests an over-conversion bug"


def test_load_fhs_matches_fhs_datalibrary_glucose():
    """``load_fhs`` and ``DataLibrary().get('FHS').load()`` agree on glucose."""
    df = load.load_fhs()
    geo = DataLibrary().get("FHS").load()

    common_ids = df.index.intersection(geo.clinical.index)
    assert len(common_ids) > 0

    # Both paths convert mg/dL to mmol/L. They should agree to floating-point
    # precision since they use the same conversion factor.
    sample_id = common_ids[0]
    assert (
        abs(
            df.loc[sample_id, "glucose"]
            - geo.clinical.loc[sample_id, "glucose"]
        )
        < 1e-6
    )


def verify_expected_columns(df):
    actual_columns = set(df.columns.to_list())
    missing_columns = set(load_columns) - actual_columns
    extra_columns = actual_columns - set(load_columns)
    assert (
        len(missing_columns) == 0
    ), f"Missing expected columns: {missing_columns} \n Found extra columns: {extra_columns}"


# Run the test
if __name__ == "__main__":
    pytest.main([__file__])
