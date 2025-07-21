import pytest
import pandas as pd
import numpy as np
from biolearn.model_gallery import ModelGallery
from biolearn.data_library import GeoData
from biolearn.util import load_test_data_file


def test_grimage_missing_sex_column():
    """Test that missing sex column gives helpful error message."""
    sample_inputs = load_test_data_file("external/DNAmTestSet.csv")
    sample_metadata = load_test_data_file("external/testset_metadata.csv")

    # Remove sex column
    metadata_no_sex = sample_metadata.drop(columns=['sex'])
    test_data = GeoData(metadata_no_sex, sample_inputs)

    gallery = ModelGallery()
    grimage_model = gallery.get('GrimAgeV2')

    with pytest.raises(ValueError) as exc_info:
        grimage_model.predict(test_data)

    error_msg = str(exc_info.value)
    assert "GrimAge requires 'sex' column" in error_msg
    assert "SexEstimation" in error_msg


def test_grimage_missing_age_column():
    """Test that missing age column gives clear error message."""
    sample_inputs = load_test_data_file("external/DNAmTestSet.csv")
    sample_metadata = load_test_data_file("external/testset_metadata.csv")

    # Remove age column
    metadata_no_age = sample_metadata.drop(columns=['age'])
    test_data = GeoData(metadata_no_age, sample_inputs)

    gallery = ModelGallery()
    grimage_model = gallery.get('GrimAgeV2')

    with pytest.raises(ValueError) as exc_info:
        grimage_model.predict(test_data)

    error_msg = str(exc_info.value)
    assert "GrimAge requires 'age' column" in error_msg


def test_grimage_nan_sex_values():
    """Test that NaN sex values give helpful error message."""
    sample_inputs = load_test_data_file("external/DNAmTestSet.csv")
    sample_metadata = load_test_data_file("external/testset_metadata.csv")

    # Set some sex values to NaN
    test_metadata = sample_metadata.copy()
    test_metadata.loc[test_metadata.index[0], 'sex'] = np.nan
    test_data = GeoData(test_metadata, sample_inputs)

    gallery = ModelGallery()
    grimage_model = gallery.get('GrimAgeV2')

    with pytest.raises(ValueError) as exc_info:
        grimage_model.predict(test_data)

    error_msg = str(exc_info.value)
    assert "cannot process samples with unknown sex" in error_msg
    assert "SexEstimation" in error_msg


def test_grimage_sample_mismatch():
    """Test that mismatched samples between methylation and metadata give clear error."""
    sample_inputs = load_test_data_file("external/DNAmTestSet.csv")
    sample_metadata = load_test_data_file("external/testset_metadata.csv")

    # Remove one sample from metadata to create mismatch
    mismatched_metadata = sample_metadata.drop(sample_metadata.index[0])
    test_data = GeoData(mismatched_metadata, sample_inputs)

    gallery = ModelGallery()
    grimage_model = gallery.get('GrimAgeV2')

    with pytest.raises(ValueError) as exc_info:
        grimage_model.predict(test_data)

    error_msg = str(exc_info.value)
    assert "Methylation data contains samples without metadata" in error_msg
    assert "Ensure all samples in methylation matrix have corresponding metadata" in error_msg


def test_grimage_reverse_sample_mismatch():
    """Test metadata samples without methylation data give clear error."""
    sample_inputs = load_test_data_file("external/DNAmTestSet.csv")
    sample_metadata = load_test_data_file("external/testset_metadata.csv")

    # Add extra sample to metadata to create reverse mismatch
    extra_metadata = sample_metadata.copy()
    extra_metadata.loc['EXTRA_SAMPLE'] = {'age': 50, 'sex': 1}
    test_data = GeoData(extra_metadata, sample_inputs)

    gallery = ModelGallery()
    grimage_model = gallery.get('GrimAgeV2')

    with pytest.raises(ValueError) as exc_info:
        grimage_model.predict(test_data)

    error_msg = str(exc_info.value)
    assert "Metadata contains samples without methylation data" in error_msg
    assert "Ensure all samples in metadata have corresponding methylation" in error_msg
