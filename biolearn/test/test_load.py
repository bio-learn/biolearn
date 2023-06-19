from biolearn import load
from biolearn.model import load_columns
import pytest
import numpy as np


def test_fhs_columns():
    df = load.load_fhs()
    verify_expected_columns(df)


def test_nhanes_columns():
    df = load.load_nhanes(2010)
    verify_expected_columns(df)


def test_can_load_nhanes_2012():
    df = load.load_nhanes(2012)


def test_can_load_dnam():
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41169/matrix/GSE41169_series_matrix.txt.gz"
    df = load.load_dnam(
        dnam_file=url, id_row=32, age_row=46, skiprows=72
    )  # nrows=1 to make it faster
    # Verify data set is of known size
    # assert df.shape == (540, 27579) need to be more general
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
