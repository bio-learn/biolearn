import pandas as pd
import numpy as np
from biolearn.data_library import GeoData


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
