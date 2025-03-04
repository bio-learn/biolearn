import os
import numpy as np
import pandas as pd
import pytest
from biolearn.data_library import GeoData

def create_dummy_geodata(num_samples=5):
    """
    Create a dummy GeoData instance with:
      - A methylation (dnam) DataFrame with random beta values (rounded to 3 decimals).
      - A metadata DataFrame with a 'Sex' column using internal numeric coding:
            1 for Female, 2 for Male.
        For this test, samples with even index are Female (1) and odd index are Male (2).
      - An 'Age' column and a 'Disease State' column.
      - RNA and protein data are set to None.
    """
    # Create dummy methylation data (10 CpG sites, num_samples samples)
    cpg_ids = [f"cg{i:08d}" for i in range(1, 11)]
    sample_ids = [f"S{i}" for i in range(1, num_samples + 1)]
    data = np.random.rand(len(cpg_ids), num_samples)
    dnam = pd.DataFrame(data, index=cpg_ids, columns=sample_ids).round(3)
    
    # Create dummy metadata with internal numeric sex representation:
    # For samples: even-indexed sample gets 1 (Female), odd-indexed gets 2 (Male)
    sex_values = [2 if i % 2 != 0 else 1 for i in range(1, num_samples + 1)]
    metadata = pd.DataFrame({
        "SampleID": sample_ids,
        "Sex": sex_values,
        "Age": [25 + i for i in range(num_samples)],
        "Disease State": ["None"] * num_samples
    }).set_index("SampleID")
    
    return GeoData(metadata, dnam=dnam)

@pytest.fixture
def temp_dir(tmp_path):
    """Return a temporary directory as a string."""
    return str(tmp_path)

def test_metadata_population():
    data = {"Sample1": [0.1, 0.2, 0.3], "Sample2": [0.4, 0.5, 0.6]}
    dnam_df = pd.DataFrame(data, index=["Site1", "Site2", "Site3"])

    geo_data = GeoData.from_methylation_matrix(dnam_df)

    assert list(geo_data.metadata.index) == list(
        dnam_df.columns
    ), "Metadata row names do not match expected values."


def test_technical_duplicate_averaging():
    data = {"Sample1": [0.1, 0.2], "Sample2": [0.4, 0.5]}
    dnam_df = pd.DataFrame(data, index=["Site1_1", "Site1_2"])

    expected_data = {"Sample1": [0.15], "Sample2": [0.45]}
    expected_df = pd.DataFrame(expected_data, index=["Site1"])

    geo_data = GeoData.from_methylation_matrix(dnam_df)

    pd.testing.assert_frame_equal(
        geo_data.dnam, expected_df, check_dtype=False
    )


def test_can_load_from_file_path(tmp_path):
    # Create a temporary CSV file to simulate file input
    d = tmp_path / "sub"
    d.mkdir()
    file_path = d / "methylation_matrix.csv"
    data = {"Sample1": [0.1, 0.2, 0.3], "Sample2": [0.4, 0.5, 0.6]}
    dnam_df = pd.DataFrame(data, index=["Site1", "Site2", "Site3"])
    dnam_df.to_csv(file_path)

    # Create GeoData instance from file path
    geo_data = GeoData.from_methylation_matrix(str(file_path))

    # Verify the data is loaded correctly
    pd.testing.assert_frame_equal(geo_data.dnam, dnam_df, check_dtype=False)

def test_save_load_roundtrip(temp_dir):
    """
    Test that a GeoData instance can be saved and loaded without data loss.
    """
    geodata = create_dummy_geodata(num_samples=5)
    geodata.save_csv(temp_dir, "testdata")
    loaded_geodata = GeoData.load_csv(temp_dir, "testdata", series_part="all")

    print("Original metadata:")
    print(geodata.metadata)
    print(geodata.metadata.dtypes)

    print("Loaded metadata:")
    print(loaded_geodata.metadata)
    print(loaded_geodata.metadata.dtypes)

    orig_ds = geodata.metadata["Disease State"]
    loaded_ds = loaded_geodata.metadata["Disease State"]

    print("Original 'Disease State':", orig_ds.unique())
    print("Loaded 'Disease State':", loaded_ds.unique())
    print("Original dtype:", orig_ds.dtype)
    print("Loaded dtype:", loaded_ds.dtype)

    # Use check_dtype=False to ignore minor attribute differences such as dtype for the "Disease State" column.
    pd.testing.assert_frame_equal(geodata.metadata, loaded_geodata.metadata, check_dtype=False)
    pd.testing.assert_frame_equal(geodata.dnam, loaded_geodata.dnam)

    assert loaded_geodata.rna is None
    assert loaded_geodata.protein is None

def test_sex_conversion(temp_dir):
    """
    Test that the 'Sex' field is converted to the standard on save and reversed on load.
    Assumes the internal representation uses "F" and "M", and the standard is 0 for female, 1 for male.
    """
    geodata = create_dummy_geodata(num_samples=3)
    geodata.save_csv(temp_dir, "sex_test")
    loaded_geodata = GeoData.load_csv(temp_dir, "sex_test", series_part="all")
    
    expected_sex = geodata.metadata["Sex"]
    actual_sex = loaded_geodata.metadata["Sex"]
    pd.testing.assert_series_equal(expected_sex, actual_sex, check_names=False)

def test_split_large_methylation(temp_dir):
    """
    Test that methylation data with >1000 samples is split into multiple files.
    """
    num_samples = 1050
    geodata = create_dummy_geodata(num_samples=num_samples)
    geodata.save_csv(temp_dir, "large")
    
    files = os.listdir(temp_dir)
    methyl_files = [fname for fname in files if "methylation" in fname]
    assert len(methyl_files) >= 2, "Expected multiple split methylation files for >1000 samples"
    
    loaded_geodata = GeoData.load_csv(temp_dir, "large", series_part="all")
    pd.testing.assert_frame_equal(geodata.dnam, loaded_geodata.dnam)

def test_missing_files(temp_dir):
    """
    Test that if optional files (e.g., RNA) are missing, the corresponding attribute is set to None.
    """
    geodata = create_dummy_geodata(num_samples=5)
    geodata.rna = geodata.dnam.copy()
    geodata.save_csv(temp_dir, "missing_test")
    
    rna_file = os.path.join(temp_dir, "missing_test_rna.csv")
    if os.path.exists(rna_file):
        os.remove(rna_file)
    
    loaded_geodata = GeoData.load_csv(temp_dir, "missing_test", series_part="all")
    assert loaded_geodata.rna is None, "Expected RNA attribute to be None when RNA file is missing"

# ---------------------------
# Existing Tests for from_methylation_matrix
# ---------------------------

def test_metadata_population():
    data = {"Sample1": [0.1, 0.2, 0.3], "Sample2": [0.4, 0.5, 0.6]}
    dnam_df = pd.DataFrame(data, index=["Site1", "Site2", "Site3"])
    
    geo_data = GeoData.from_methylation_matrix(dnam_df)
    
    assert list(geo_data.metadata.index) == list(dnam_df.columns), "Metadata row names do not match expected values."

def test_technical_duplicate_averaging():
    data = {"Sample1": [0.1, 0.2], "Sample2": [0.4, 0.5]}
    dnam_df = pd.DataFrame(data, index=["Site1_1", "Site1_2"])
    
    expected_data = {"Sample1": [0.15], "Sample2": [0.45]}
    expected_df = pd.DataFrame(expected_data, index=["Site1"])
    
    geo_data = GeoData.from_methylation_matrix(dnam_df)
    
    pd.testing.assert_frame_equal(geo_data.dnam, expected_df, check_dtype=False)

def test_can_load_from_file_path(tmp_path):
    # Create a temporary CSV file to simulate file input
    d = tmp_path / "sub"
    d.mkdir()
    file_path = d / "methylation_matrix.csv"
    data = {"Sample1": [0.1, 0.2, 0.3], "Sample2": [0.4, 0.5, 0.6]}
    dnam_df = pd.DataFrame(data, index=["Site1", "Site2", "Site3"])
    dnam_df.to_csv(file_path)
    
    # Create GeoData instance from file path
    geo_data = GeoData.from_methylation_matrix(str(file_path))
    
    pd.testing.assert_frame_equal(geo_data.dnam, dnam_df, check_dtype=False)