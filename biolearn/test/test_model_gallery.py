import pytest
from biolearn.util import get_data_file
from biolearn.model import LinearMethylationModel
from biolearn.model_gallery import ModelGallery


def no_transform(_):
    return _


def test_methylation_model_can_report_sites():
    coefficient_file = get_data_file("Horvath1.csv")
    test_model = LinearMethylationModel(coefficient_file, no_transform)

    some_expected_cpgs = ["cg00945507", "cg04528819", "cg14727952"]
    expected_cpg_count = 353

    result_sites = test_model.methylation_sites()
    assert expected_cpg_count == len(
        result_sites
    ), f"Expected count {expected_cpg_count} but got {len(result_sites)}"
    assert all(
        cpg in result_sites for cpg in some_expected_cpgs
    ), f"Some expected CpGs not found in result"


# Sample Model Data for testing
sample_data = {
    "Horvathv1": {
        "year": 2013,
        "species": "Human",
        "tissue": "Multi-tissue",
        "source": "https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115",
        "model": {
            "type": "LinearMethylationModel",
            "file": "Horvath1.csv",
            "transform": lambda sum: sum + 0.696,
        },
    },
    "SampleModel2": {
        "year": 2020,
        "species": "Human",
        "tissue": "Brain",
        "source": "https://samplelink.com",
        "model": {
            "type": "LinearMethylationModel",
            "file": "Horvath2.csv",
            "transform": lambda sum: sum + 1.0,
        },
    },
}


def test_init_with_bad_data():
    # This example assumes a bad data structure, adjust accordingly
    bad_data = {"ModelName": "SomeStringInsteadOfDict"}
    with pytest.raises(
        ValueError, match="Expected dictionary for model definition, got"
    ):
        ModelGallery(bad_data)


def test_get_model_by_name():
    gallery = ModelGallery(sample_data)
    model = gallery.get("Horvathv1")
    assert model.metadata["year"] == 2013


def test_get_nonexistent_model_by_name():
    gallery = ModelGallery(sample_data)
    with pytest.raises(KeyError, match="Model not found:"):
        gallery.get("NonExistentModel")


def test_search_with_no_parameters():
    gallery = ModelGallery(sample_data)
    results = gallery.search()
    assert len(results) == len(sample_data)
    assert "Horvathv1" in results
    assert "SampleModel2" in results


def test_search_by_parameters():
    gallery = ModelGallery(sample_data)
    results = gallery.search(species="Human", tissue="Brain")
    assert len(results) == 1
    assert "SampleModel2" in results
