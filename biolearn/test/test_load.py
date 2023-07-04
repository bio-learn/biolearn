from biolearn import load
from biolearn.model import load_columns
import pytest
import numpy as np
import os


def test_fhs_columns():
    df = load.load_fhs()
    verify_expected_columns(df)


def test_nhanes_columns():
    df = load.load_nhanes(2010)
    verify_expected_columns(df)


def test_can_load_nhanes_2012():
    df = load.load_nhanes(2012)


def test_can_load_dnam():
    script_dir = os.path.dirname(
        __file__
    )  # get the directory of the current script
    data_file_path = os.path.join(
        script_dir, "data", "geo_dnam_test_file"
    )  # build the path to the data file
    df = load.load_dnam(
        dnam_file=data_file_path, id_row=32, age_row=46, skiprows=72)
    # Verify data set is of known size
    assert df.shape == (5, 38)
    assert "age" in df.columns.to_list()
    assert all(np.issubdtype(df[col].dtype, np.number) for col in df.columns)



def verify_expected_columns(df):
    actual_columns = set(df.columns.to_list())
    missing_columns = set(load_columns) - actual_columns
    extra_columns = actual_columns - set(load_columns)
    assert (
        len(missing_columns) == 0
    ), f"Missing expected columns: {missing_columns} \n Found extra columns: {extra_columns}"


# Run the test
if __name__ == "__main__":
    pytest.main([__file__])
