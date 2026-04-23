import pytest
from biolearn import model
from biolearn.model_gallery import ModelGallery

gallery = ModelGallery()


@pytest.mark.parametrize(
    "model_name, model_entry", model.model_definitions.items()
)
def test_required_features_interface(model_name, model_entry):
    """Every model must implement required_features() with the correct shape."""
    model_type = model_entry["model"]["type"]

    # Skip types that can't be instantiated without special setup
    if model_type in ["NotImplemented", "HurdleAPIModel"]:
        pytest.skip(f"Model type {model_type} requires special setup")

    model_class = getattr(model, model_type)
    instance = model_class.from_definition(model_entry)

    result = instance.required_features()

    # Validate return format
    assert isinstance(
        result, dict
    ), f"{model_name}: required_features() must return a dict"
    assert (
        "layer" in result
    ), f"{model_name}: required_features() must include 'layer'"
    assert (
        "features" in result
    ), f"{model_name}: required_features() must include 'features'"
    assert (
        "metadata" in result
    ), f"{model_name}: required_features() must include 'metadata'"

    assert isinstance(
        result["layer"], str
    ), f"{model_name}: 'layer' must be a string"
    assert isinstance(
        result["features"], list
    ), f"{model_name}: 'features' must be a list"
    assert isinstance(
        result["metadata"], list
    ), f"{model_name}: 'metadata' must be a list"

    # Layer should be one of the known GeoData layers
    valid_layers = {
        "dnam",
        "rna",
        "protein_alamar",
        "protein_olink",
        "clinical",
    }
    assert (
        result["layer"] in valid_layers
    ), f"{model_name}: 'layer' must be one of {valid_layers}, got '{result['layer']}'"


@pytest.mark.parametrize(
    "model_name, model_entry", model.model_definitions.items()
)
def test_required_features_consistency_with_methylation_sites(
    model_name, model_entry
):
    """For dnam models, required_features() features should match methylation_sites()."""
    model_type = model_entry["model"]["type"]

    if model_type in ["NotImplemented", "HurdleAPIModel"]:
        pytest.skip(f"Model type {model_type} requires special setup")

    model_class = getattr(model, model_type)
    instance = model_class.from_definition(model_entry)

    if not hasattr(instance, "methylation_sites"):
        pytest.skip(f"{model_name} does not have methylation_sites()")

    result = instance.required_features()
    sites = instance.methylation_sites()

    if result["layer"] == "dnam" and sites:
        assert set(result["features"]) == set(
            sites
        ), f"{model_name}: required_features() features should match methylation_sites()"
