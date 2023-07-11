import pytest
import os
from math import isclose
import pandas as pd
from biolearn import clock


def load_test_data_file(relative_path):
    script_dir = os.path.dirname(
        __file__
    )  # get the directory of the current script
    data_file_path = os.path.join(
        script_dir, "data", relative_path
    )  # build the path to the data file
    test_sample = pd.read_csv(data_file_path, index_col=0)
    return test_sample


def load_test_sample():
    return load_test_data_file("dnam-test-sample.csv")


def test_horvathv1():
    sample_data = load_test_sample()
    expected = 62.85
    actual = clock.single_sample_clock(
        clock.horvathv1, sample_data.transpose()
    )
    assert isclose(actual, expected, abs_tol=1e-2)


def test_hannum():
    sample_data = load_test_sample()
    # This seems very strange
    expected = -3.05
    actual = clock.single_sample_clock(
        clock.hannum, sample_data.transpose()
    )
    assert isclose(actual, expected, abs_tol=1e-2)


def test_phenoage():
    sample_data = load_test_sample()
    expected = 61.09
    actual = clock.single_sample_clock(
        clock.phenoage, sample_data.transpose()
    )
    assert isclose(actual, expected, abs_tol=1e-2)


# Run the test
if __name__ == "__main__":
    pytest.main([__file__])
