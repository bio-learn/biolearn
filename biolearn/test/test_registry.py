import pandas as pd
import pytest
from biolearn.clinical.registry import BIOMARKER_REGISTRY
from biolearn.clinical.convert import convert_units, validate_ranges


class TestBiomarkerRegistry:
    def test_known_biomarkers_not_empty(self):
        assert len(BIOMARKER_REGISTRY) > 0

    def test_get_albumin(self):
        entry = BIOMARKER_REGISTRY.get("albumin")
        # Canonical unit matches what NHANES exposes via LBDSALSI
        assert entry["unit"] == "g/L"
        assert "range" in entry
        assert "conversions" in entry

    def test_get_unknown_raises(self):
        with pytest.raises(KeyError, match="Unknown biomarker"):
            BIOMARKER_REGISTRY.get("nonexistent_biomarker")

    def test_canonical_unit(self):
        # Both glucose and creatinine canonical units match NHANES SI columns
        assert BIOMARKER_REGISTRY.canonical_unit("glucose") == "mmol/L"
        assert BIOMARKER_REGISTRY.canonical_unit("creatinine") == "umol/L"

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
        assert names == sorted(names)

    def test_source_preset_nhanes(self):
        preset = BIOMARKER_REGISTRY.get_source_preset("nhanes")
        assert isinstance(preset, dict)
        # NHANES values are pre-canonicalized in load_nhanes, so no conversions
        assert preset == {}

    def test_source_preset_fhs(self):
        preset = BIOMARKER_REGISTRY.get_source_preset("fhs")
        # FHS Period 1 ships glucose in mg/dL; canonical is mmol/L
        assert preset.get("glucose") == "mg/dL"

    def test_source_preset_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown source preset"):
            BIOMARKER_REGISTRY.get_source_preset("fake_source")

    def test_ukbiobank_preset_is_not_registered(self):
        """We removed the UK Biobank preset because we have not validated it."""
        with pytest.raises(ValueError, match="Unknown source preset"):
            BIOMARKER_REGISTRY.get_source_preset("ukbiobank")


class TestConvertUnits:
    def test_glucose_mg_per_dl_to_mmol_per_l(self):
        df = pd.DataFrame({"glucose": [90.0]}, index=["P1"])
        result = convert_units(df, units={"glucose": "mg/dL"})
        assert abs(result.loc["P1", "glucose"] - 90.0 * 0.05551) < 1e-6

    def test_hdl_cholesterol_mg_per_dl_to_mmol_per_l(self):
        df = pd.DataFrame({"hdl_cholesterol": [50.0]}, index=["P1"])
        result = convert_units(df, units={"hdl_cholesterol": "mg/dL"})
        assert abs(result.loc["P1", "hdl_cholesterol"] - 50.0 / 38.67) < 1e-3

    def test_creatinine_mg_per_dl_to_umol_per_l(self):
        # Canonical is umol/L; the conversion from mg/dL multiplies by 88.42
        df = pd.DataFrame({"creatinine": [1.0]}, index=["P1"])
        result = convert_units(df, units={"creatinine": "mg/dL"})
        assert abs(result.loc["P1", "creatinine"] - 88.42) < 1e-3

    def test_albumin_g_per_dl_to_g_per_l(self):
        # Canonical is g/L; the conversion from g/dL multiplies by 10
        df = pd.DataFrame({"albumin": [4.2]}, index=["P1"])
        result = convert_units(df, units={"albumin": "g/dL"})
        assert abs(result.loc["P1", "albumin"] - 42.0) < 1e-6

    def test_no_conversion_returns_copy(self):
        df = pd.DataFrame({"albumin": [4.2]}, index=["P1"])
        result = convert_units(df)
        pd.testing.assert_frame_equal(result, df)
        assert result is not df  # should be a different object

    def test_already_canonical_no_change(self):
        # Canonical creatinine is umol/L; passing umol/L should be a no-op
        df = pd.DataFrame({"creatinine": [60.0]}, index=["P1"])
        result = convert_units(df, units={"creatinine": "umol/L"})
        assert result.loc["P1", "creatinine"] == 60.0

    def test_unknown_unit_raises(self):
        df = pd.DataFrame({"creatinine": [1.0]}, index=["P1"])
        with pytest.raises(ValueError, match="No conversion"):
            convert_units(df, units={"creatinine": "fake_unit"})

    def test_missing_column_skipped(self):
        df = pd.DataFrame({"albumin": [42.0]}, index=["P1"])
        result = convert_units(df, units={"creatinine": "mg/dL"})
        assert "albumin" in result.columns

    def test_fhs_source_preset(self):
        df = pd.DataFrame({"glucose": [90.0]}, index=["P1"])
        result = convert_units(df, source_units="fhs")
        assert abs(result.loc["P1", "glucose"] - 90.0 * 0.05551) < 1e-6


class TestValidateRanges:
    def test_in_range_no_warnings(self):
        # Albumin canonical unit is g/L; 40 is mid-range
        df = pd.DataFrame({"albumin": [40.0]}, index=["P1"])
        result = validate_ranges(df, warn=False)
        assert len(result) == 0

    def test_out_of_range_detected(self):
        # 0.1 g/L is well below the (10, 60) g/L range
        df = pd.DataFrame({"albumin": [0.1]}, index=["P1"])
        result = validate_ranges(df, warn=False)
        assert "albumin" in result
        assert result["albumin"] == 1

    def test_unknown_columns_ignored(self):
        df = pd.DataFrame({"unknown_col": [999]}, index=["P1"])
        result = validate_ranges(df, warn=False)
        assert len(result) == 0
