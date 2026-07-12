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
