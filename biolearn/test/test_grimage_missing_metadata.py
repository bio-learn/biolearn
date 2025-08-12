import pytest
import pandas as pd
import numpy as np
from biolearn.model_gallery import ModelGallery
from biolearn.data_library import GeoData
from biolearn.util import get_test_data_file


def get_test_data():
    return GeoData.load_csv(get_test_data_file("testset/"), "testset")


def test_grimage_missing_sex_column():
    """Test that missing sex column gives helpful error message."""
    test_data = get_test_data()
    test_data.metadata = test_data.metadata.drop(columns=["sex"])

    gallery = ModelGallery()
    grimage_model = gallery.get("GrimAgeV2")

    with pytest.raises(ValueError) as exc_info:
        grimage_model.predict(test_data)

    error_msg = str(exc_info.value)
    assert "GrimAge requires 'sex' column" in error_msg
    assert "SexEstimation" in error_msg


def test_grimage_missing_age_column():
    """Test that missing age column gives clear error message."""
    test_data = get_test_data()
    test_data.metadata = test_data.metadata.drop(columns=["age"])

    gallery = ModelGallery()
    grimage_model = gallery.get("GrimAgeV2")

    with pytest.raises(ValueError) as exc_info:
        grimage_model.predict(test_data)

    error_msg = str(exc_info.value)
    assert "GrimAge requires 'age' column" in error_msg


def test_grimage_nan_sex_values():
    """Test that NaN sex values give helpful error message."""
    test_data = get_test_data()
    test_data.metadata.loc[test_data.metadata.index[0], "sex"] = np.nan

    gallery = ModelGallery()
    grimage_model = gallery.get("GrimAgeV2")

    with pytest.raises(ValueError) as exc_info:
        grimage_model.predict(test_data)

    error_msg = str(exc_info.value)
    assert "cannot process samples with unknown sex" in error_msg
    assert "SexEstimation" in error_msg


def test_grimage_sample_mismatch():
    """Test that mismatched samples between methylation and metadata give clear error."""
    test_data = get_test_data()
    # Remove one sample from metadata to create mismatch
    test_data.metadata = test_data.metadata.iloc[1:]

    gallery = ModelGallery()
    grimage_model = gallery.get("GrimAgeV2")

    with pytest.raises(ValueError) as exc_info:
        grimage_model.predict(test_data)

    error_msg = str(exc_info.value)
    assert "Methylation data contains samples without metadata" in error_msg
    assert (
        "Ensure all samples in methylation matrix have corresponding metadata"
        in error_msg
    )
