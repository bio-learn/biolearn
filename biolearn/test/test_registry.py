import pandas as pd
import pytest
from biolearn.clinical.registry import BIOMARKER_REGISTRY
from biolearn.clinical.convert import convert_units, validate_ranges


class TestBiomarkerRegistry:
    def test_known_biomarkers_not_empty(self):
        assert len(BIOMARKER_REGISTRY) > 0

    def test_get_albumin(self):
        entry = BIOMARKER_REGISTRY.get("albumin")
        assert entry["unit"] == "g/dL"
        assert "range" in entry
        assert "conversions" in entry

    def test_get_unknown_raises(self):
        with pytest.raises(KeyError, match="Unknown biomarker"):
            BIOMARKER_REGISTRY.get("nonexistent_biomarker")

    def test_canonical_unit(self):
        assert BIOMARKER_REGISTRY.canonical_unit("glucose") == "mmol/L"
        assert BIOMARKER_REGISTRY.canonical_unit("creatinine") == "mg/dL"

    def test_valid_range(self):
        lo, hi = BIOMARKER_REGISTRY.valid_range("albumin")
        assert lo < hi

    def test_contains(self):
        assert "albumin" in BIOMARKER_REGISTRY
        assert "fake_marker" not in BIOMARKER_REGISTRY

    def test_known_biomarkers_list(self):
        names = BIOMARKER_REGISTRY.known_biomarkers()
        assert isinstance(names, list)
        assert "albumin" in names
        assert names == sorted(names)  # should be sorted

    def test_source_preset_nhanes(self):
        preset = BIOMARKER_REGISTRY.get_source_preset("nhanes")
        assert isinstance(preset, dict)

    def test_source_preset_ukbiobank(self):
        preset = BIOMARKER_REGISTRY.get_source_preset("ukbiobank")
        assert "creatinine" in preset
        assert preset["creatinine"] == "umol/L"

    def test_source_preset_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown source preset"):
            BIOMARKER_REGISTRY.get_source_preset("fake_source")


class TestConvertUnits:
    def test_creatinine_umol_to_mg(self):
        df = pd.DataFrame({"creatinine": [88.42]}, index=["P1"])
        result = convert_units(df, units={"creatinine": "umol/L"})
        assert abs(result.loc["P1", "creatinine"] - 1.0) < 0.01

    def test_albumin_g_per_l_to_g_per_dl(self):
        df = pd.DataFrame({"albumin": [42.0]}, index=["P1"])
        result = convert_units(df, units={"albumin": "g/L"})
        assert abs(result.loc["P1", "albumin"] - 4.2) < 0.01

    def test_no_conversion_returns_copy(self):
        df = pd.DataFrame({"albumin": [4.2]}, index=["P1"])
        result = convert_units(df)
        pd.testing.assert_frame_equal(result, df)
        assert result is not df  # should be a different object

    def test_already_canonical_no_change(self):
        df = pd.DataFrame({"creatinine": [1.0]}, index=["P1"])
        result = convert_units(df, units={"creatinine": "mg/dL"})
        assert result.loc["P1", "creatinine"] == 1.0

    def test_unknown_unit_raises(self):
        df = pd.DataFrame({"creatinine": [1.0]}, index=["P1"])
        with pytest.raises(ValueError, match="No conversion"):
            convert_units(df, units={"creatinine": "fake_unit"})

    def test_missing_column_skipped(self):
        df = pd.DataFrame({"albumin": [4.2]}, index=["P1"])
        result = convert_units(df, units={"creatinine": "umol/L"})
        assert "albumin" in result.columns

    def test_source_preset(self):
        df = pd.DataFrame(
            {"creatinine": [88.42], "albumin": [42.0]}, index=["P1"]
        )
        result = convert_units(df, source_units="ukbiobank")
        assert abs(result.loc["P1", "creatinine"] - 1.0) < 0.01
        assert abs(result.loc["P1", "albumin"] - 4.2) < 0.01


class TestValidateRanges:
    def test_in_range_no_warnings(self):
        df = pd.DataFrame({"albumin": [4.0]}, index=["P1"])
        result = validate_ranges(df, warn=False)
        assert len(result) == 0

    def test_out_of_range_detected(self):
        df = pd.DataFrame({"albumin": [0.1]}, index=["P1"])  # below range
        result = validate_ranges(df, warn=False)
        assert "albumin" in result
        assert result["albumin"] == 1

    def test_unknown_columns_ignored(self):
        df = pd.DataFrame({"unknown_col": [999]}, index=["P1"])
        result = validate_ranges(df, warn=False)
        assert len(result) == 0
