import pytest

from biolearn.features import (
    RequiredFeatures,
    MissingFeaturesError,
    VALID_LAYERS,
)


def test_required_features_behaves_as_object():
    rf = RequiredFeatures("dnam", ("cg1", "cg2"))
    assert rf.layer == "dnam"
    assert rf.features == ("cg1", "cg2")
    assert rf.metadata == ()


def test_required_features_behaves_as_mapping():
    rf = RequiredFeatures("dnam", ("cg1", "cg2"))
    assert rf["layer"] == "dnam"
    assert rf["features"] == ["cg1", "cg2"]
    assert rf["metadata"] == []
    assert set(rf.keys()) == {"layer", "features", "metadata"}
    assert dict(rf) == {
        "layer": "dnam",
        "features": ["cg1", "cg2"],
        "metadata": [],
    }


def test_metadata_key_always_present_when_empty():
    rf = RequiredFeatures("clinical", ("albumin",))
    assert "metadata" in rf
    assert rf["metadata"] == []


def test_missing_features_error_is_value_error():
    assert issubclass(MissingFeaturesError, ValueError)


def test_valid_layers_contains_expected():
    assert "dnam" in VALID_LAYERS
    assert "clinical" in VALID_LAYERS


import types

import pandas as pd

from biolearn.features import validate_required_features


class _StubModel:
    def __init__(self, req, name=None):
        self._req = req
        self.details = {"name": name} if name else {}

    def required_features(self):
        return self._req


def _geo(dnam=None, rna=None, clinical=None, metadata=None):
    return types.SimpleNamespace(
        dnam=dnam,
        rna=rna,
        protein_olink=None,
        protein_alamar=None,
        clinical=clinical,
        metadata=metadata,
    )


def test_validate_passes_when_dnam_present():
    dnam = pd.DataFrame({"S1": [0.1, 0.2]}, index=["cg1", "cg2"])
    model = _StubModel(RequiredFeatures("dnam", ("cg1", "cg2")))
    assert validate_required_features(model, _geo(dnam=dnam)) is None


def test_validate_dnam_missing_raises_with_cpg_message():
    dnam = pd.DataFrame({"S1": [0.1]}, index=["cg1"])
    model = _StubModel(RequiredFeatures("dnam", ("cg1", "cg2")), name="Clock")
    with pytest.raises(
        MissingFeaturesError, match=r"Missing required CpG sites.*cg2"
    ):
        validate_required_features(model, _geo(dnam=dnam))


def test_validate_clinical_uses_columns_axis():
    # clinical is samples-as-rows, biomarkers-as-columns
    clinical = pd.DataFrame({"albumin": [4.0], "glucose": [5.0]}, index=["S1"])
    model = _StubModel(RequiredFeatures("clinical", ("albumin", "crp")))
    with pytest.raises(MissingFeaturesError, match=r"crp"):
        validate_required_features(model, _geo(clinical=clinical))


def test_validate_metadata_missing_raises():
    dnam = pd.DataFrame({"S1": [0.1]}, index=["cg1"])
    meta = pd.DataFrame({"sex": ["m"]}, index=["S1"])
    model = _StubModel(RequiredFeatures("dnam", ("cg1",), ("age",)))
    with pytest.raises(
        MissingFeaturesError, match=r"Missing required metadata.*age"
    ):
        validate_required_features(model, _geo(dnam=dnam, metadata=meta))


def test_validate_missing_layer_frame_is_empty_set():
    # layer present but geo_data has None for it -> all features missing
    model = _StubModel(RequiredFeatures("dnam", ("cg1",)))
    with pytest.raises(MissingFeaturesError):
        validate_required_features(model, _geo(dnam=None))


from biolearn.model import (
    LinearMethylationModel,
    LinearTranscriptomicModel,
)


def _linear_methylation_model():
    coeffs = pd.DataFrame(
        {"CoefficientTraining": [0.5, 0.5, 1.0]},
        index=["cg1", "cg2", "intercept"],
    )
    return LinearMethylationModel(coeffs, transform=lambda x: x, name="TestLM")


def test_linear_methylation_required_features():
    model = _linear_methylation_model()
    rf = model.required_features()
    assert rf.layer == "dnam"
    assert set(rf.features) == {"cg1", "cg2"}
    assert rf.metadata == ()
    assert "intercept" not in rf.features


def test_linear_transcriptomic_layer_is_rna():
    coeffs = pd.DataFrame(
        {"CoefficientTraining": [0.5, 1.0]}, index=["GENE1", "intercept"]
    )
    model = LinearTranscriptomicModel(coeffs, transform=lambda x: x)
    rf = model.required_features()
    assert rf.layer == "rna"
    assert set(rf.features) == {"GENE1"}


def test_linear_methylation_predict_raises_on_missing_cpg():
    model = _linear_methylation_model()
    dnam = pd.DataFrame({"S1": [0.1]}, index=["cg1"])
    geo = _geo(dnam=dnam)
    with pytest.raises(
        MissingFeaturesError, match=r"Missing required CpG sites.*cg2"
    ):
        model.predict(geo)


def test_linear_methylation_missing_cpg_is_value_error():
    # Backward compatibility: existing `except ValueError` still catches it.
    model = _linear_methylation_model()
    dnam = pd.DataFrame({"S1": [0.1]}, index=["cg1"])
    with pytest.raises(ValueError):
        model.predict(_geo(dnam=dnam))
