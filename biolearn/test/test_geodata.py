import os
import numpy as np
import pandas as pd
import pytest
from biolearn.data_library import GeoData
from biolearn.util import get_test_data_file


def create_dummy_geodata(num_samples=5):
    """
    Create a dummy GeoData instance with:
      - A methylation (dnam) DataFrame with random beta values (rounded to 3 decimals).
      - A metadata DataFrame with a 'Sex' column using standard numeric coding:
            0 for Female, 1 for Male.
        For this test, samples with even index are Female (0) and odd index are Male (1).
      - An 'Age' column and a 'Disease State' column.
      - RNA and protein data are set to None.
    """
    # Create dummy methylation data (10 CpG sites, num_samples samples)
    cpg_ids = [f"cg{i:08d}" for i in range(1, 11)]
    sample_ids = [f"S{i}" for i in range(1, num_samples + 1)]
    data = np.random.rand(len(cpg_ids), num_samples)
    dnam = pd.DataFrame(data, index=cpg_ids, columns=sample_ids).round(3)

    sex_values = [1 if i % 2 != 0 else 0 for i in range(1, num_samples + 1)]
    metadata = pd.DataFrame(
        {
            "SampleID": sample_ids,
            "Sex": sex_values,
            "Age": [25 + i for i in range(num_samples)],
            "Disease State": ["None"] * num_samples,
        }
    ).set_index("SampleID")

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
    pd.testing.assert_frame_equal(geo_data.dnam, dnam_df, check_dtype=False)


def test_save_load_roundtrip(temp_dir):
    """
    Test that a GeoData instance can be saved and loaded without data loss.
    """
    geodata = create_dummy_geodata(num_samples=5)
    geodata.save_csv(temp_dir, "testdata")
    loaded_geodata = GeoData.load_csv(temp_dir, "testdata", series_part="all")

    pd.testing.assert_frame_equal(
        geodata.metadata, loaded_geodata.metadata, check_dtype=False
    )
    pd.testing.assert_frame_equal(geodata.dnam, loaded_geodata.dnam)
    assert loaded_geodata.rna is None
    assert loaded_geodata.protein_alamar is None


def test_sex_conversion(temp_dir):
    """
    Test that the 'Sex' field is converted to the standard on save and reversed on load.
    """
    geodata = create_dummy_geodata(num_samples=3)
    geodata.save_csv(temp_dir, "sex_test")
    loaded_geodata = GeoData.load_csv(temp_dir, "sex_test", series_part="all")
    expected_sex = geodata.metadata["Sex"]
    actual_sex = loaded_geodata.metadata["Sex"]
    pd.testing.assert_series_equal(expected_sex, actual_sex, check_names=False)


def test_sex_conversion_on_disk(temp_dir):
    """
    Test that the metadata file on disk conforms to the spec:
    The 'Sex' field should be saved in standard format: 0 for female and 1 for male.
    Given our internal coding (1 for Female, 2 for Male), the expected values on disk for 3 samples are:
      - For sample S1 (internal 2) -> 1
      - For sample S2 (internal 1) -> 0
      - For sample S3 (internal 2) -> 1
    """
    geodata = create_dummy_geodata(num_samples=3)
    geodata.save_csv(temp_dir, "sex_disk")
    metadata_file = os.path.join(temp_dir, "sex_disk_metadata.csv")
    saved_metadata = pd.read_csv(
        metadata_file, index_col=0, keep_default_na=False
    )
    expected_sex_on_disk = [
        1,
        0,
        1,
    ]  # Internal 2 becomes 1, internal 1 becomes 0
    actual_sex_on_disk = saved_metadata["Sex"].tolist()
    assert (
        actual_sex_on_disk == expected_sex_on_disk
    ), f"Expected Sex on disk {expected_sex_on_disk} but got {actual_sex_on_disk}"


def test_split_large_methylation(temp_dir):
    """
    Test that methylation data with >1000 samples is split into multiple files.
    """
    num_samples = 1050
    geodata = create_dummy_geodata(num_samples=num_samples)
    geodata.save_csv(temp_dir, "large")
    files = os.listdir(temp_dir)
    methyl_files = [fname for fname in files if "methylation" in fname]
    assert (
        len(methyl_files) >= 2
    ), "Expected multiple split methylation files for >1000 samples"
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
    loaded_geodata = GeoData.load_csv(
        temp_dir, "missing_test", series_part="all"
    )
    assert (
        loaded_geodata.rna is None
    ), "Expected RNA attribute to be None when RNA file is missing"


def test_unspecified_proteomic_file(temp_dir):
    """
    Test that if an unspecified proteomic file is present, the function raises
    an exception with the phrase "unspecified source proteomic file found".
    """
    geodata = create_dummy_geodata(num_samples=5)
    geodata.save_csv(temp_dir, "test_data")

    # Create a proteomic file in the directory
    proteomic_file = os.path.join(temp_dir, "test_data_protein.csv")
    with open(proteomic_file, "w") as f:
        f.write("dummy,proteomic,data\n")
        f.write("value1,value2,value3\n")

    # Assert that loading raises an exception with the expected phrase
    with pytest.raises(Exception) as exc_info:
        loaded_geodata = GeoData.load_csv(
            temp_dir, "test_data", series_part="all"
        )

    assert (
        "unspecified source proteomic file found"
        in str(exc_info.value).lower()
    ), f"Expected exception message to contain 'unspecified source proteomic file found', but got: {exc_info.value}"


def test_optional_field_content(temp_dir):
    """
    Test that when RNA and protein data are provided, they are correctly saved and loaded.
    """
    geodata = create_dummy_geodata(num_samples=5)
    geodata.rna = geodata.dnam.copy()
    geodata.protein_alamar = geodata.dnam.copy() * 10
    geodata.save_csv(temp_dir, "optional")
    loaded_geodata = GeoData.load_csv(temp_dir, "optional", series_part="all")
    pd.testing.assert_frame_equal(geodata.rna, loaded_geodata.rna)
    pd.testing.assert_frame_equal(
        geodata.protein_alamar, loaded_geodata.protein_alamar
    )


def test_missing_metadata_file(temp_dir):
    """
    Test that if the metadata file is missing, the loaded GeoData.metadata is None.
    """
    geodata = create_dummy_geodata(num_samples=5)
    geodata.save_csv(temp_dir, "missing_meta")
    metadata_file = os.path.join(temp_dir, "missing_meta_metadata.csv")
    if os.path.exists(metadata_file):
        os.remove(metadata_file)
    loaded_geodata = GeoData.load_csv(
        temp_dir, "missing_meta", series_part="all"
    )
    assert (
        loaded_geodata.metadata is None
    ), "Expected metadata to be None when metadata file is missing"


def test_missing_methylation_file(temp_dir):
    """
    Test that if the methylation file(s) are missing, the loaded GeoData.dnam is None.
    """
    geodata = create_dummy_geodata(num_samples=5)
    geodata.save_csv(temp_dir, "missing_methyl")
    for fname in os.listdir(temp_dir):
        if fname.startswith("missing_methyl_methylation"):
            os.remove(os.path.join(temp_dir, fname))
    loaded_geodata = GeoData.load_csv(
        temp_dir, "missing_methyl", series_part="all"
    )
    assert (
        loaded_geodata.dnam is None
    ), "Expected dnam to be None when methylation file is missing"


def test_invalid_series_part_parameter(temp_dir):
    """
    Test that an invalid series_part parameter (non-numeric and not 'all') raises a ValueError.
    """
    geodata = create_dummy_geodata(num_samples=5)
    geodata.save_csv(temp_dir, "invalid_series")
    with pytest.raises(ValueError):
        GeoData.load_csv(temp_dir, "invalid_series", series_part="invalid")


def test_boundary_for_splitting_methylation(temp_dir):
    """
    Test the boundary condition for splitting methylation data.
    For exactly 1000 samples, only one file should be created.
    For 1001 samples, multiple files should be created.
    """
    # Exactly 1000 samples.
    geodata_1000 = create_dummy_geodata(num_samples=1000)
    geodata_1000.save_csv(temp_dir, "boundary_1000")
    files_1000 = os.listdir(temp_dir)
    methyl_files_1000 = [
        f for f in files_1000 if f.startswith("boundary_1000_methylation")
    ]
    assert (
        len(methyl_files_1000) == 1
    ), "Expected one methylation file for exactly 1000 samples"

    # Clear out for next test.
    for f in methyl_files_1000:
        os.remove(os.path.join(temp_dir, f))

    # 1001 samples.
    geodata_1001 = create_dummy_geodata(num_samples=1001)
    geodata_1001.save_csv(temp_dir, "boundary_1001")
    files_1001 = os.listdir(temp_dir)
    methyl_files_1001 = [
        f for f in files_1001 if f.startswith("boundary_1001_methylation")
    ]
    assert (
        len(methyl_files_1001) >= 2
    ), "Expected multiple methylation files for 1001 samples"


def test_loading_specific_methylation_part(temp_dir):
    """
    Test that loading a specific methylation part (e.g., part 2) returns the correct subset of columns.
    """
    num_samples = 1500  # This should split into at least 2 parts.
    geodata = create_dummy_geodata(num_samples=num_samples)
    geodata.save_csv(temp_dir, "specific_part")

    # Determine how many parts were created.
    all_files = os.listdir(temp_dir)
    part_files = sorted(
        [f for f in all_files if f.startswith("specific_part_methylation")]
    )
    assert (
        len(part_files) >= 2
    ), "Expected at least 2 methylation files for 1500 samples"

    # Load part 2 specifically.
    loaded_part2 = GeoData.load_csv(temp_dir, "specific_part", series_part=2)
    start = 1000
    end = num_samples  # Assuming part 2 covers columns 1001 to 1500.
    expected_dnam = geodata.dnam.iloc[:, start:end]
    pd.testing.assert_frame_equal(expected_dnam, loaded_part2.dnam)


def test_file_overwriting_behavior(temp_dir):
    """
    Test that saving a GeoData instance twice (with modifications in between)
    overwrites the previous files and the loaded data reflects the changes.
    """
    geodata = create_dummy_geodata(num_samples=5)
    geodata.save_csv(temp_dir, "overwrite")

    # Modify metadata: change Age for all samples.
    geodata.metadata["Age"] = geodata.metadata["Age"] + 10
    geodata.save_csv(temp_dir, "overwrite")

    loaded_geodata = GeoData.load_csv(temp_dir, "overwrite", series_part="all")
    pd.testing.assert_frame_equal(
        geodata.metadata, loaded_geodata.metadata, check_dtype=False
    )
    pd.testing.assert_frame_equal(geodata.dnam, loaded_geodata.dnam)


def test_type_preservation_in_metadata(temp_dir):
    """
    Test that numeric metadata fields such as 'Age' remain numeric (int or float) after save and load.
    """
    geodata = create_dummy_geodata(num_samples=5)
    geodata.save_csv(temp_dir, "type_preservation")
    loaded_geodata = GeoData.load_csv(
        temp_dir, "type_preservation", series_part="all"
    )
    assert pd.api.types.is_numeric_dtype(
        loaded_geodata.metadata["Age"]
    ), "Expected Age column to be numeric"


def test_dnam_columns_numeric_and_single_nan():
    file_path = get_test_data_file("")
    data = GeoData.load_csv(file_path, "example")

    # Check all columns in dnam are numeric
    for col in data.dnam.columns:
        assert pd.api.types.is_numeric_dtype(
            data.dnam[col]
        ), f"Column {col} is not numeric."

    # Ensure there is exactly one NaN value in the dnam DataFrame
    nan_count = data.dnam.isna().sum().sum()
    assert nan_count == 1, f"Expected 1 NaN value, but found {nan_count}."
